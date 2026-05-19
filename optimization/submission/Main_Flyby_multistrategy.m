clear; clc

%% =====================================================================
%  MAIN_FLYBY_MULTISTRATEGY  --  SCORE-HUNT VERSION
%  ---------------------------------------------------------------------
%  Purpose
%  -------
%  This script runs a multi-strategy beam search for a propellant-free
%  gravity-assist tour in the Mythara system. The trajectory is modeled as
%  a sequence of Lambert arcs joined by zero-duration patched-conic flybys.
%  Each node in the search tree represents one partial tour: body sequence,
%  absolute flyby epochs, Lambert branch flags, incoming V-infinity at the
%  current body, raw science score, official score, and a heuristic rank.
%
%  The main goal of this HEAVY/HARD version is not a quick feasibility run.
%  It is intended to explore many possible body sequences and transfer-time
%  combinations in order to find high-scoring tours, especially tours that
%  include repeated valuable flybys of Ashkara, Varuun, and Xerion.
%
%  Dynamical model
%  ---------------
%  1. Planet states are generated with KM_solver.m from bodies.mat.
%  2. Transfers between two consecutive flybys are solved with
%     lambert_universal.m.
%  3. At each intermediate body, the incoming and outgoing hyperbolic excess
%     velocities must have the same magnitude.
%  4. The required turning angle is checked against the patched-conic flyby
%     model, using the allowed periapsis-altitude range [0.1, 100] body radii.
%  5. Nyxar is massless. If it is allowed as an interior flyby, the full
%     V-infinity vector must be continuous, not only its magnitude.
%  6. Each Lambert arc is also checked against the minimum heliocentric
%     perihelion constraint.
%
%  Search strategy
%  ---------------
%  The search is a deterministic/stochastic hybrid beam search:
%
%     Stage 1: generate all valid two-body seed arcs over a first-TOF grid.
%     Stage 2: expand the seeds once before applying tight beam pruning.
%     Stage 3: repeatedly expand each node by trying every target body, both
%              Lambert branches, and several optimized time-of-flight roots.
%     Stage 4: rank nodes using current score, score density, remaining time,
%              future potential, body diversity, and penalties for slow tours.
%     Stage 5: prune with bucket protection so that the beam does not collapse
%              into only one family of similar sequences.
%     Stage 6: save a global archive, export the best CSV files, and optionally
%              run a local time-refinement pass on elite nodes.
%
%  Scoring
%  -------
%  The script keeps two scores for every node:
%
%     J_raw      : science score without the Grand Tour multiplier.
%     J_official : science score with optional Grand Tour multiplier.
%
%  The active objective is selected per configuration with:
%
%     cfg.USE_GRAND_TOUR_BONUS = true/false
%
%  In most raw-score searches this is false, because the priority is to
%  maximize the score validated by the web site without relying on the bonus.
%
%  How to use
%  ----------
%  Put this file in the optimization folder and run:
%
%       Main_Flyby_multistrategy
%
%  Required files on the MATLAB path or nearby folders:
%
%       KM_solver.m
%       lambert_universal.m
%       stumpff.m
%       solution_to_csv_patched.m
%       bodies.mat
%
%  Output
%  ------
%  The script creates results_multistrategy/ containing:
%
%       archive_global_final.mat   : full archive and best node
%       archive_summary.csv        : readable summary of all retained nodes
%       csv/                       : exported CSV candidates for validation
%       run_*.mat                  : result of each individual configuration
%
%  Notes for tuning
%  ----------------
%  This is a heavy search. Runtime can be many hours or days. For debugging,
%  reduce BEAM_WIDTH, SEED_BEAM_WIDTH, MAX_DEPTH, TOF_GRID_SHORT_N,
%  TOF_GRID_LONG_N, N_TOP_INTERVALS, and N_TOP_TOF.
%
%  This documented version includes the v3 bug fixes:
%     - safe candidate struct creation in refine_next_tof_multi();
%     - row-vector-safe index concatenation in prune_frontier_diverse().
%
%  =====================================================================

% Locate this script, infer the project root, and add all required folders
% to the MATLAB path. This avoids accidental use of shadowed functions.
script_dir = fileparts(mfilename('fullpath'));
[root, paths] = setup_paths(script_dir);

bodies_file = find_bodies_file(root, paths, script_dir);
load(bodies_file,'bodies')

mu = bodies.mu_star;
ys = bodies.year_sec;
if isfield(bodies,'AU')
    AU = bodies.AU;
else
    AU = 149597870.691;
end

fprintf('\nUsing:\n')
which KM_solver
which lambert_universal
which stumpff
which solution_to_csv_patched
fprintf('bodies.mat = %s\n', bodies_file)

% Build global numerical/search options and then build the list of run
% configurations. Each configuration changes the body ordering, beam width,
% TOF grid, stochastic mode, and Grand Tour option.
opts = default_opts(root, paths, bodies_file);

% Start/resize the MATLAB parallel pool before the heavy search.
% The expensive parts of this script are parallelized with parfor over
% independent seed arcs, parent-node expansions, and local refinements.
opts = start_parallel_pool(opts);

cfgs = build_configs(opts);

if ~isfolder(opts.OUTDIR)
    mkdir(opts.OUTDIR)
end

% Global archive accumulates good nodes from every configuration.
% best_global is always selected by J_raw, even when a configuration uses
% the Grand Tour bonus internally for ranking.
archive_global = [];
best_global = [];

fprintf('\n=== MAIN MULTISTRATEGY SEARCH ===\n')
fprintf('Output dir: %s\n', opts.OUTDIR)
fprintf('Number of configs: %d\n', numel(cfgs))

for icfg = 1:numel(cfgs)
    cfg = cfgs(icfg);
    rng(cfg.RNG_SEED);

    fprintf('\n============================================================\n')
    fprintf('RUN %02d/%02d | %s\n', icfg, numel(cfgs), cfg.NAME)
    fprintf('============================================================\n')
    fprintf('BODY_IDS = [%s]\n', sprintf('%d ', cfg.BODY_IDS))
    fprintf('MAX_DEPTH=%d | BEAM=%d | STOCH=%d | GT_BONUS=%d\n', ...
        cfg.MAX_DEPTH, cfg.BEAM_WIDTH, cfg.USE_STOCHASTIC_BEAM, cfg.USE_GRAND_TOUR_BONUS)

    res = run_one_search(cfg, opts, bodies, mu, ys, AU);

    archive_global = merge_archive(archive_global, res.archive, opts.MAX_GLOBAL_ARCHIVE);

    if isempty(best_global) || res.best.J_raw > best_global.J_raw
        best_global = res.best;
    end

    save_run_result(res, cfg, opts);

    if cfg.N_REFINE_TOP > 0 && ~isempty(res.archive)
        fprintf('\n--- local refinement on top %d nodes ---\n', cfg.N_REFINE_TOP)
        elites = take_top_unique(res.archive, cfg.N_REFINE_TOP);
        refined_cells = cell(1, numel(elites));
        if opts.USE_PARALLEL && numel(elites) > 1
            parfor k = 1:numel(elites)
                refined_cells{k} = refine_sequence_times(elites(k), cfg, opts, bodies, mu, ys, AU);
            end
        else
            for k = 1:numel(elites)
                refined_cells{k} = refine_sequence_times(elites(k), cfg, opts, bodies, mu, ys, AU);
            end
        end

        refined = concat_node_cells(refined_cells);
        if ~isempty(refined)
            [~, ir] = max([refined.J_raw]);
            if isempty(best_global) || refined(ir).J_raw > best_global.J_raw
                best_global = refined(ir);
                fprintf('REFINE NEW BEST raw=%.9f t=%.3f seq=[%s]\n', ...
                    best_global.J_raw, best_global.t(end), sprintf('%d ', best_global.seq))
            end
        end
        archive_global = merge_archive(archive_global, refined, opts.MAX_GLOBAL_ARCHIVE);
    end
end

archive_global = prune_archive_unique(archive_global, opts.MAX_GLOBAL_ARCHIVE);
archive_global = sort_nodes(archive_global, 'raw');

summary_file = fullfile(opts.OUTDIR, 'archive_summary.csv');
write_archive_summary(archive_global, summary_file, bodies)

fprintf('\n=== GLOBAL BEST RAW ===\n')
disp_best(best_global)
print_sequence_details(best_global, bodies, mu, ys)

export_top_csvs(archive_global, opts, opts.N_EXPORT_CSV)

save(fullfile(opts.OUTDIR, 'archive_global_final.mat'), ...
    'archive_global', 'best_global', 'cfgs', 'opts')

fprintf('\nDONE.\n')
fprintf('Best raw = %.9f\n', best_global.J_raw)
fprintf('Summary  = %s\n', summary_file)
fprintf('Folder   = %s\n', opts.OUTDIR)

%% =====================================================================
function opts = start_parallel_pool(opts)
%START_PARALLEL_POOL Start a MATLAB parallel pool using the requested cores.
% This script is configured here for a 6-core machine. Each worker is forced to one
% computational thread to avoid CPU oversubscription.

if ~isfield(opts, 'USE_PARALLEL') || ~opts.USE_PARALLEL
    fprintf('Parallel execution disabled. Running serially.\n')
    return
end

if isempty(ver('parallel'))
    warning('Parallel Computing Toolbox not found. Running serially.')
    opts.USE_PARALLEL = false;
    return
end

nWorkers = max(1, opts.N_WORKERS);
try
    c = parcluster('local');
    nWorkers = min(nWorkers, c.NumWorkers);
catch
    % Keep requested value if parcluster is unavailable for some reason.
end

pool = gcp('nocreate');
if isempty(pool)
    fprintf('Starting parallel pool with %d workers...\n', nWorkers)
    pool = parpool('local', nWorkers);
elseif pool.NumWorkers ~= nWorkers
    fprintf('Restarting parallel pool with %d workers...\n', nWorkers)
    delete(pool)
    pool = parpool('local', nWorkers);
else
    fprintf('Using existing parallel pool with %d workers.\n', pool.NumWorkers)
end

try
    maxNumCompThreads(opts.THREADS_PER_WORKER);
    pctRunOnAll('maxNumCompThreads(1);');
catch ME
    warning('Could not force one thread per worker: %s', ME.message)
end

try
    pctRunOnAll('rehash toolboxcache;');
catch
end

fprintf('Parallel mode active: %d workers, %d thread per worker.\n', ...
    pool.NumWorkers, opts.THREADS_PER_WORKER)

end

%% =====================================================================
function opts = default_opts(root, paths, bodies_file)
%DEFAULT_OPTS Define global constants and all user-tunable controls.
%
% =====================================================================
%                         USER TUNING PANEL
% =====================================================================
% This is the main section to modify when you want to increase the score.
% The most important knobs are:
%
%   1) BEAM_WIDTH / SEED_BEAM_WIDTH in build_configs()
%      -> larger values keep more candidate trajectories alive.
%
%   2) TOF_GRID_SHORT_N / TOF_GRID_LONG_N / N_TOP_INTERVALS / N_TOP_TOF
%      -> larger values test more transfer-time alternatives.
%
%   3) RAW_PROTECT_FRAC
%      -> protects a fraction of the beam using J_raw directly, so good
%         high-score families are not removed only because the heuristic
%         rank is temporarily lower.
%
%   4) Ranking weights below
%      -> tune how much the search prefers score density, heavy bodies,
%         diversity, future potential, short transfers, etc.
%
% Suggested tuning strategy:
%   - First increase BEAM_WIDTH and RAW_PROTECT_FRAC.
%   - Then increase TOF grid and N_TOP_TOF.
%   - Only after that modify rank weights.
% =====================================================================

opts.ROOT = root;
opts.PATHS = paths;
opts.BODIES_FILE = bodies_file;
opts.OUTDIR = fullfile(root, 'results_multistrategy');

%% =====================================================================
%  A) PARALLEL EXECUTION SETTINGS
% =====================================================================
% Keep this active for your 6-core machine.
% THREADS_PER_WORKER = 1 avoids oversubscription: N workers x 1 thread is
% usually better than workers where each one tries to use many threads.

opts.USE_PARALLEL = true;
opts.N_WORKERS = min(6, feature('numcores'));
opts.THREADS_PER_WORKER = 1;

%% =====================================================================
%  B) GLOBAL PHYSICAL / MISSION CONSTRAINTS
% =====================================================================
% Be careful with these: they define the validity of the trajectory search.
% You normally do NOT tune these just to increase score.

opts.MAX_TIME = 300;                % maximum mission duration [yr]
opts.MAX_FLYBYS_PER_BODY = 13;      % website/scoring limit per body

opts.DIRS = [0 1];                  % Lambert branches to test
opts.TOF_MIN = 1.5;                 % minimum time of flight [yr]
opts.TOF_MAX = 160;                 % maximum time of flight [yr]
opts.VINFTOL = 1e-3;                % strict V-infinity closure tolerance [km/s]
opts.SOFT_VINFTOL_FACTOR = 3.0;     % softer tolerance during candidate generation

opts.MIN_H_OVER_R = 0.1;            % minimum flyby altitude / body radius
opts.MAX_H_OVER_R = 100.0;          % maximum flyby altitude / body radius
opts.MIN_PERIHELION_AU = 0.01;      % heliocentric perihelion safety constraint
opts.DEFAULT_PERIHELION_AU = 0.05;

%% =====================================================================
%  C) TRANSFER-TIME SEARCH SETTINGS  <--- IMPORTANT TO TUNE
% =====================================================================
% These parameters control how many TOF roots are searched for each possible
% next leg. Increasing them improves exploration but increases runtime.
%
% Light values:       36 / 34 / 5 / 3
% Recommended strong: 50 / 46 / 8 / 4
% Very strong:        60 / 55 / 10 / 5
% Extreme:            70 / 65 / 12 / 6

opts.TOF_GRID_SHORT_N = 52;         % grid points from TOF_MIN to TOF_SHORT_CUT_YR
opts.TOF_GRID_LONG_N  = 46;         % grid points from TOF_SHORT_CUT_YR to TOF_MAX
opts.TOF_SHORT_CUT_YR = 35;         % split between dense short grid and long grid
opts.N_TOP_INTERVALS = 7;            % number of promising TOF intervals refined by fminbnd
opts.N_TOP_TOF = 4;                  % number of accepted TOF candidates kept per leg

% Soft cost terms used only while searching for TOF candidates.
% COST_TURN_W penalizes impossible flyby turn angles during TOF search.
% COST_TOF_W mildly penalizes long transfers during TOF search.
% If you increase COST_TOF_W too much, the algorithm may over-prefer short
% arcs and miss high-score long transfers.
opts.COST_TURN_W = 1e-4;
opts.COST_TOF_W  = 2e-5;

%% =====================================================================
%  D) PRUNING / ARCHIVE SETTINGS  <--- VERY IMPORTANT TO TUNE
% =====================================================================
% RAW_PROTECT_FRAC is one of the most useful changes for your case.
% It forces pruning to keep a fraction of nodes selected by J_raw directly,
% not only by the heuristic rank.
%
% Recommended:
%   0.15 = conservative
%   0.25 = strong / good default
%   0.35 = aggressive, keeps many high raw-score nodes

opts.RAW_PROTECT_FRAC = 0.30;
opts.RAW_PROTECT_MIN = 80;

% Bucket protection preserves trajectory diversity.
% Smaller TIME_BIN_YR and VINF_BIN create more buckets -> more diversity.
% Larger KEEP_PER_BUCKET keeps more variants from each family.
opts.TIME_BIN_YR = 10;
opts.VINF_BIN = 1.5;
opts.KEEP_PER_BUCKET = 3;

% Larger archive keeps more globally good solutions for later export/refine.
opts.MAX_GLOBAL_ARCHIVE = 60000;
opts.N_EXPORT_CSV = 40;

%% =====================================================================
%  E) STOCHASTIC BEAM SETTINGS
% =====================================================================
% Used only in configs with USE_STOCHASTIC_BEAM = true.
% ELITE_FRAC = deterministic top part of the beam.
% SOFTMAX_T controls randomness: higher = more random, lower = greedier.

opts.ELITE_FRAC = 0.30;
opts.SOFTMAX_T = 18.0;

%% =====================================================================
%  F) RANKING PARAMETERS / HEURISTIC WEIGHTS  <--- WHERE TO TUNE WEIGHTS
% =====================================================================
% rank_node_smart() uses these values. The final objective is still J_raw,
% but rank decides which partial trajectories survive beam pruning.
%
% If you want to push score higher, usually increase:
%   RANK_WEIGHT_RAW_SCORE        -> more greedy on current score
%   RANK_WEIGHT_SCORE_RATE       -> favors score per year
%   RANK_FUTURE_SCALE            -> favors paths with future potential
%   RANK_WEIGHT_HEAVY_UNIQUE     -> favors visiting high-value bodies
%
% If the search gets stuck in repeated low-value loops, increase penalties:
%   PEN_LOW_4_W, PEN_LOW_5_W, PEN_LAG_W, PEN_LATE_W

opts.TARGET_MEAN_TOF = 18;
opts.TARGET_RECENT_TOF = 15;
opts.REPEAT_FRAC = 1/3;
opts.MIN_INSERT_GAP = 5.0;

% Body groups used by the heuristic.
% HEAVY_IDS are bodies the ranking should protect/push.
% LOWVALUE_IDS are bodies that should not dominate the beam too much.
opts.HEAVY_IDS = [10 9 8 7];
opts.LOWVALUE_IDS = [4 5 1000];

% Future-potential estimate.
opts.RANK_F_REF = 0.72;
opts.RANK_FUTURE_SCALE = 0.24;

% Main positive rank weights.
opts.RANK_WEIGHT_RAW_SCORE = 1.0;       % multiplier of current J_raw
opts.RANK_WEIGHT_SCORE_RATE = 18.0;     % multiplier of J_raw / elapsed time
opts.RANK_WEIGHT_HEAVY_UNIQUE = 2.5;    % reward for unique heavy bodies visited
opts.RANK_WEIGHT_ALL_UNIQUE = 0.8;      % reward for total body diversity

% Penalty weights.
opts.PEN_SLOW_W = 2.0;                  % penalty for average TOF above target
opts.PEN_RECENT_SLOW_W = 1.5;           % penalty for recent TOF above target
opts.PEN_LOW_4_W = 3;                 % penalty after too many visits to body 4
opts.PEN_LOW_5_W = 1.5;                 % penalty after too many visits to body 5
opts.PEN_NYXAR_W = 5.0;                 % penalty for Nyxar in raw-score runs
opts.PEN_LATE_W = 30.0;                 % penalty for consuming too much total time
opts.PEN_LAG_W = 4.5;                   % penalty if number of flybys is low vs time

%% =====================================================================
%  G) LOCAL REFINEMENT SETTINGS
% =====================================================================
% Local refinement keeps the body sequence fixed and perturbs the leg times.
% It can improve final candidates, but it is not the main global search.

opts.REFINE_RESTARTS = 3;
opts.REFINE_SWEEPS = 2;
opts.REFINE_STEP_FRAC = 0.04;
opts.REFINE_MIN_STEP_YR = 0.03;

end


%% =====================================================================
function cfgs = build_configs(opts)
%BUILD_CONFIGS Create the set of independent search strategies.
%
% =====================================================================
%                      CONFIGURATION TUNING PANEL
% =====================================================================
% Each make_cfg line defines one complete search run.
%
% make_cfg(base, name, BODY_IDS, use_GT_bonus, stochastic, first_TOF_grid,
%          BEAM_WIDTH, MAX_DEPTH, SEED_BEAM_WIDTH, RNG_SEED, N_REFINE_TOP)
%
% What to tune here:
%
%   BODY_IDS
%       Order of bodies tested during expansion. It is not a physical
%       constraint, but it affects tie-breaking and pruning.
%
%   BEAM_WIDTH
%       Number of partial trajectories kept after each depth.
%       Bigger = better exploration but slower.
%       Good values: 1500 quick, 3000 strong, 4500 very strong, 6000 extreme.
%
%   SEED_BEAM_WIDTH
%       Number of initial two-body seeds kept before depth-3 expansion.
%       Usually 3x-5x BEAM_WIDTH.
%
%   MAX_DEPTH
%       Maximum number of bodies in the sequence. Do not increase too much:
%       the 300-year limit is usually more restrictive than depth.
%
%   N_REFINE_TOP
%       Number of elite sequences locally refined after each config.
% =====================================================================

base = base_cfg(opts);
cfgs = base([]);

cfgs(end+1) = make_cfg(base, 'exploit_249_spine_det', ...
    [9 8 6 5 7 4], false, false, ...
    [1.5:0.25:30 31:0.75:85], ...
    1200, 24, 5000, 249, 0);
cfgs(end).N_REFINE_TOP = 0;
cfgs(end).USE_ENTRY_FILTER = false;
cfgs(end).ALLOW_NYXAR_INTERIOR_PASS = false;
cfgs(end).ALLOW_NYXAR_TERMINAL_ONLY = false;
cfgs(end).USE_GRAND_TOUR_BONUS = false;
cfgs(end).GT_BONUS = false;
% Medium run variant (commented): BEAM=1600, SEED_BEAM=7000, MAX_DEPTH=28
% Heavy  run variant (commented): BEAM=2200, SEED_BEAM=10000, MAX_DEPTH=30

cfgs(end+1) = make_cfg(base, 'exploit_249_spine_stoch', ...
    [9 8 6 5 7 4], false, true, ...
    [1.5:0.25:30 31:0.75:85], ...
    1400, 24, 6000, 250, 0);
cfgs(end).N_REFINE_TOP = 0;
cfgs(end).USE_ENTRY_FILTER = false;
cfgs(end).ALLOW_NYXAR_INTERIOR_PASS = false;
cfgs(end).ALLOW_NYXAR_TERMINAL_ONLY = false;
cfgs(end).USE_GRAND_TOUR_BONUS = false;
cfgs(end).GT_BONUS = false;
% Medium run variant (commented): BEAM=1600, SEED_BEAM=7000, MAX_DEPTH=28
% Heavy  run variant (commented): BEAM=2200, SEED_BEAM=10000, MAX_DEPTH=30

cfgs(end+1) = make_cfg(base, 'exploit_241_249_hybrid', ...
    [9 8 6 7 5 4], false, true, ...
    [1.5:0.20:28 29:0.50:80], ...
    1400, 25, 6000, 251, 0);
cfgs(end).N_REFINE_TOP = 0;
cfgs(end).USE_ENTRY_FILTER = false;
cfgs(end).ALLOW_NYXAR_INTERIOR_PASS = false;
cfgs(end).ALLOW_NYXAR_TERMINAL_ONLY = false;
cfgs(end).USE_GRAND_TOUR_BONUS = false;
cfgs(end).GT_BONUS = false;
% Medium run variant (commented): BEAM=1600, SEED_BEAM=7000, MAX_DEPTH=28
% Heavy  run variant (commented): BEAM=2200, SEED_BEAM=10000, MAX_DEPTH=30

end


%% =====================================================================
function cfg = base_cfg(opts)
%BASE_CFG Default settings inherited by every individual configuration.
% A configuration created by make_cfg() overrides the fields that matter for
% that specific run.

cfg.NAME = 'base';
cfg.BODY_IDS = [10 9 8 7 6 5 4];
cfg.USE_GRAND_TOUR_BONUS = false;
cfg.GT_BONUS = false;
cfg.USE_STOCHASTIC_BEAM = false;
cfg.USE_ENTRY_FILTER = false;
cfg.ALLOW_NYXAR_INTERIOR_PASS = false;
cfg.ALLOW_NYXAR_TERMINAL_ONLY = false;

cfg.MAX_TIME = opts.MAX_TIME;
cfg.MAX_DEPTH = 36;
cfg.BEAM_WIDTH = 900;
cfg.SEED_BEAM_WIDTH = 4000;
cfg.SEED_PRUNE_DEPTH = 3;
cfg.TOF_FIRST_GRID = [2:0.5:80];
cfg.DIRS = opts.DIRS;
cfg.RNG_SEED = 1;
cfg.N_REFINE_TOP = 0;

end

%% =====================================================================
function cfg = make_cfg(base, name, body_ids, use_gt, stochastic, tof_grid, beam_width, max_depth, seed_beam, rng_seed, n_refine)
%MAKE_CFG Helper constructor for one search configuration.

cfg = base;
cfg.NAME = name;
cfg.BODY_IDS = body_ids;
cfg.USE_GRAND_TOUR_BONUS = use_gt;
cfg.GT_BONUS = use_gt;
cfg.USE_STOCHASTIC_BEAM = stochastic;
cfg.TOF_FIRST_GRID = tof_grid;
cfg.BEAM_WIDTH = beam_width;
cfg.MAX_DEPTH = max_depth;
cfg.SEED_BEAM_WIDTH = seed_beam;
cfg.RNG_SEED = rng_seed;
cfg.N_REFINE_TOP = n_refine;

end

%% =====================================================================
function res = run_one_search(cfg, opts, bodies, mu, ys, AU)
%RUN_ONE_SEARCH Execute one complete beam-search configuration.
% The function starts from all valid pair seeds, expands them to depth 3
% before tight pruning, then continues layer by layer until MAX_DEPTH or until
% no valid children remain.

frontier = generate_pair_seeds(cfg, opts, bodies, mu, ys, AU);
spine_frontier = build_spine_seed_nodes_249(cfg, opts, bodies, mu, ys, AU);
fprintf('Injected spine seeds: %d\n', numel(spine_frontier))
frontier = merge_archive(frontier, spine_frontier, cfg.SEED_BEAM_WIDTH);

if isempty(frontier)
    warning('No initial seed found for %s.', cfg.NAME)
    res.best = [];
    res.archive = [];
    return
end

fprintf('Initial seed count before prune: %d\n', numel(frontier))
frontier = prune_frontier_diverse(frontier, cfg.SEED_BEAM_WIDTH, cfg, opts);
print_family249_diag(frontier)
fprintf('Initial seed count after prune:  %d\n', numel(frontier))
print_best(frontier, 5)

archive = frontier;
[~, ib] = max([frontier.J_raw]);
best = frontier(ib);
current_depth = 2;

% Delayed pruning: expand seeds once before main pruning.
if cfg.SEED_PRUNE_DEPTH >= 3
    fprintf('\n=== SEED EXPANSION TO DEPTH 3 ===\n')
    children = expand_frontier(frontier, cfg, opts, bodies, mu, ys, AU);
    if ~isempty(children)
        archive = merge_archive(archive, children, opts.MAX_GLOBAL_ARCHIVE);
        [~, ibc] = max([children.J_raw]);
        if children(ibc).J_raw > best.J_raw
            best = children(ibc);
        end
        frontier = prune_frontier_diverse(children, cfg.BEAM_WIDTH, cfg, opts);
        print_family249_diag(frontier)
        current_depth = 3;
        print_best(frontier, 8)
    else
        fprintf('No depth-3 children; continue from seeds.\n')
    end
end

for depth = current_depth+1:cfg.MAX_DEPTH
    fprintf('\n=== DEPTH %d ===\n', depth)

    children = expand_frontier(frontier, cfg, opts, bodies, mu, ys, AU);

    if isempty(children)
        fprintf('No children found. Stop run.\n')
        break
    end

    archive = merge_archive(archive, children, opts.MAX_GLOBAL_ARCHIVE);

    [~, ibc] = max([children.J_raw]);
    if children(ibc).J_raw > best.J_raw
        best = children(ibc);
        fprintf('NEW BEST raw=%.9f t=%.3f seq=[%s]\n', ...
            best.J_raw, best.t(end), sprintf('%d ', best.seq))
    end

    frontier = prune_frontier_diverse(children, cfg.BEAM_WIDTH, cfg, opts);
    print_family249_diag(frontier)
    print_best(frontier, 8)
end

res.best = best;
res.archive = prune_archive_unique(archive, opts.MAX_GLOBAL_ARCHIVE);

fprintf('\n--- RUN BEST | %s ---\n', cfg.NAME)
disp_best(best)

end

%% =====================================================================
function frontier = generate_pair_seeds(cfg, opts, bodies, mu, ys, AU)
%GENERATE_PAIR_SEEDS Build the initial frontier from all two-body arcs.
% Parallel version: every (departure body, arrival body, first TOF, Lambert
% branch) seed is independent, so the seed generation is distributed over
% the active parallel pool when opts.USE_PARALLEL is true.

body_ids = cfg.BODY_IDS(:).';
tof_grid = cfg.TOF_FIRST_GRID(:).';
dirs = cfg.DIRS(:).';

tasks = zeros(0,4);
for a = body_ids
    for b = body_ids
        if b == a
            continue
        end
        for tof = tof_grid
            for dir = dirs
                tasks(end+1,:) = [a b tof dir]; %#ok<AGROW>
            end
        end
    end
end

nTasks = size(tasks,1);
node_cells = cell(1,nTasks);
fail_msgs = cell(1,nTasks);

if opts.USE_PARALLEL && nTasks > 1
    parfor it = 1:nTasks
        [node_cells{it}, fail_msgs{it}] = generate_one_seed(tasks(it,:), cfg, opts, bodies, mu, ys, AU);
    end
else
    for it = 1:nTasks
        [node_cells{it}, fail_msgs{it}] = generate_one_seed(tasks(it,:), cfg, opts, bodies, mu, ys, AU);
    end
end

frontier = concat_node_cells(node_cells);

% Print only a few failures to keep the command window readable.
fail_msgs = fail_msgs(~cellfun(@isempty, fail_msgs));
max_fail_print = min(20, numel(fail_msgs));
for k = 1:max_fail_print
    fprintf('%s\n', fail_msgs{k})
end
if numel(fail_msgs) > max_fail_print
    fprintf('... %d additional seed failures suppressed.\n', numel(fail_msgs) - max_fail_print)
end

end

%% =====================================================================
function spine_nodes = build_spine_seed_nodes_249(cfg, opts, bodies, mu, ys, AU)
%BUILD_SPINE_SEED_NODES_249 Inject fixed-time prefixes from validated families.

spine_nodes = [];

seq249 = [9 8 6 5 7 5 6 5 7 8 6 5 7 4 5 6 7 4 5 6 7 5 6 7 5 8 9];
seq241 = [9 8 6 7 5 6 5 7 8 6 5 7 5 4 7 6 5 4 6 7 5 4 5 7 5 6 8 9];

t241 = [0.000000 21.500000 36.260887 46.788639 54.207124 ...
        74.427167 81.582433 86.030145 98.778374 110.487324 ...
        116.005565 128.765155 147.659562 152.208761 166.223785 ...
        168.916678 177.073665 179.224759 191.084367 193.581698 ...
        211.013380 217.331253 232.780923 242.848582 257.932975 ...
        260.609410 270.130614 291.812898];

lens249 = [8 12 16 20 24 27];
lens241 = [8 12 16 20 24 28];

[t249, dirs249] = find_sequence_data_in_mat(seq249, opts);
[~, dirs241] = find_sequence_data_in_mat(seq241, opts);

if isempty(t249)
    fprintf('TODO build_spine_seed_nodes_249: t249 not found in MAT archive; injecting only seq241 prefixes.\n')
end

for L = lens241
    seq_pref = seq241(1:L);
    t_pref = t241(1:L);
    dirs_hint = [];
    if numel(dirs241) >= L-1
        dirs_hint = dirs241(1:L-1);
    end
    node_pref = build_prefix_node_from_spine(seq_pref, t_pref, dirs_hint, cfg, opts, bodies, mu, ys, AU);
    if ~isempty(node_pref)
        spine_nodes = append_node(spine_nodes, node_pref);
    end
end

if ~isempty(t249) && numel(t249) >= max(lens249)
    for L = lens249
        seq_pref = seq249(1:L);
        t_pref = t249(1:L);
        dirs_hint = [];
        if numel(dirs249) >= L-1
            dirs_hint = dirs249(1:L-1);
        end
        node_pref = build_prefix_node_from_spine(seq_pref, t_pref, dirs_hint, cfg, opts, bodies, mu, ys, AU);
        if ~isempty(node_pref)
            spine_nodes = append_node(spine_nodes, node_pref);
        end
    end
end

spine_nodes = prune_archive_unique(spine_nodes, max(1, cfg.SEED_BEAM_WIDTH));

end

%% =====================================================================
function [t_seq, dirs_seq] = find_sequence_data_in_mat(seq_target, opts)

t_seq = [];
dirs_seq = [];

mat_files = {};
if isfield(opts, 'OUTDIR')
    mat_files{end+1} = fullfile(opts.OUTDIR, 'archive_global_final.mat'); %#ok<AGROW>
    run_mats = dir(fullfile(opts.OUTDIR, 'run_*.mat'));
    for k = 1:numel(run_mats)
        mat_files{end+1} = fullfile(run_mats(k).folder, run_mats(k).name); %#ok<AGROW>
    end
end

if isfield(opts, 'ROOT') && isfolder(opts.ROOT)
    arch_mats = dir(fullfile(opts.ROOT, '**', 'archive_global_final.mat'));
    for k = 1:numel(arch_mats)
        mat_files{end+1} = fullfile(arch_mats(k).folder, arch_mats(k).name); %#ok<AGROW>
    end

    run_mats_root = dir(fullfile(opts.ROOT, '**', 'run_*.mat'));
    for k = 1:numel(run_mats_root)
        mat_files{end+1} = fullfile(run_mats_root(k).folder, run_mats_root(k).name); %#ok<AGROW>
    end
end

mat_files = unique(mat_files, 'stable');

for im = 1:numel(mat_files)
    f = mat_files{im};
    if ~isfile(f)
        continue
    end

    try
        S = load(f);
    catch
        continue
    end

    nodes = collect_nodes_from_loaded_mat(S);
    for i = 1:numel(nodes)
        if ~isfield(nodes(i), 'seq') || ~isfield(nodes(i), 't')
            continue
        end
        if isequal(nodes(i).seq, seq_target) && numel(nodes(i).t) == numel(seq_target)
            t_seq = nodes(i).t;
            if isfield(nodes(i), 'dirs')
                dirs_seq = nodes(i).dirs;
            end
            return
        end
    end
end

end

%% =====================================================================
function nodes = collect_nodes_from_loaded_mat(S)

nodes = [];

if isfield(S, 'archive_global') && ~isempty(S.archive_global)
    nodes = [nodes S.archive_global]; %#ok<AGROW>
end
if isfield(S, 'best_global') && ~isempty(S.best_global)
    nodes = [nodes S.best_global]; %#ok<AGROW>
end

if isfield(S, 'res') && isstruct(S.res)
    if isfield(S.res, 'archive') && ~isempty(S.res.archive)
        nodes = [nodes S.res.archive]; %#ok<AGROW>
    end
    if isfield(S.res, 'best') && ~isempty(S.res.best)
        nodes = [nodes S.res.best]; %#ok<AGROW>
    end
end

end

%% =====================================================================
function node = build_prefix_node_from_spine(seq_pref, t_pref, dirs_hint, cfg, opts, bodies, mu, ys, AU)
%BUILD_PREFIX_NODE_FROM_SPINE Rebuild a fixed prefix with validated Lambert dirs.

node = [];
n = numel(seq_pref);

if n < 2 || numel(t_pref) ~= n
    return
end

if abs(t_pref(1)) > 1e-10 || any(diff(t_pref) <= 0) || t_pref(end) > cfg.MAX_TIME
    return
end

if ~isempty(dirs_hint) && numel(dirs_hint) >= n-1
    node = try_build_node_from_times(seq_pref, t_pref, dirs_hint(1:n-1), cfg, opts, bodies, mu, ys, AU);
    if ~isempty(node)
        return
    end
end

dirs_pool = unique([cfg.DIRS(:).' 0 1], 'stable');
partial = [];

for dir_first = dirs_pool
    try
        [vdep, varr, ~, vp2, r1, ~] = leg_vel(seq_pref(1), seq_pref(2), t_pref(1), t_pref(2), dir_first, bodies, mu, ys);
        rp_AU = perihelion_AU(r1, vdep, mu, AU);
        if rp_AU < opts.MIN_PERIHELION_AU
            continue
        end
        arr_vinf_last = varr(:) - vp2(:);
        nd = pack_node(seq_pref(1:2), t_pref(1:2), dir_first, arr_vinf_last, NaN, NaN, rp_AU);
        nd = update_node_scores_and_rank(nd, cfg, opts, bodies, mu, ys);
        partial = append_node(partial, nd);
    catch
        continue
    end
end

if isempty(partial)
    return
end

prefix_beam = 128;

for leg = 2:n-1
    target = seq_pref(leg+1);
    t1 = t_pref(leg);
    t2 = t_pref(leg+1);
    next_nodes = [];

    for ip = 1:numel(partial)
        parent = partial(ip);

        if parent.seq(end) ~= seq_pref(leg)
            continue
        end
        if abs(parent.t(end) - t1) > 1e-8
            continue
        end

        for dir_leg = dirs_pool
            try
                [vdep, varr, vp1, vp2, r1, ~] = leg_vel(parent.seq(end), target, t1, t2, dir_leg, bodies, mu, ys);
                rp_AU = perihelion_AU(r1, vdep, mu, AU);
                if rp_AU < opts.MIN_PERIHELION_AU
                    continue
                end

                vinf_in = parent.arr_vinf_last(:);
                vinf_out = vdep(:) - vp1(:);
                [okfly, hR, delta_deg] = flyby_ok(parent.seq(end), vinf_in, vinf_out, bodies, opts, cfg);
                if ~okfly
                    continue
                end

                arr_vinf_last = varr(:) - vp2(:);
                child = pack_node([parent.seq target], [parent.t t2], [parent.dirs dir_leg], ...
                    arr_vinf_last, hR, delta_deg, rp_AU);
                child = update_node_scores_and_rank(child, cfg, opts, bodies, mu, ys);
                next_nodes = append_node(next_nodes, child);
            catch
                continue
            end
        end
    end

    if isempty(next_nodes)
        return
    end

    next_nodes = prune_archive_unique(next_nodes, prefix_beam);
    next_nodes = sort_nodes(next_nodes, 'rank');
    partial = next_nodes(1:min(prefix_beam, numel(next_nodes)));
end

if isempty(partial)
    return
end

partial = sort_nodes(partial, 'raw');
node = partial(1);

end

%% =====================================================================
function [node, fail_msg] = generate_one_seed(task, cfg, opts, bodies, mu, ys, AU)
%GENERATE_ONE_SEED Worker-safe helper for one initial two-body arc.

node = [];
fail_msg = '';

a = task(1);
b = task(2);
tof = task(3);
dir = task(4);

t1 = 0;
t2 = tof;

if t2 > cfg.MAX_TIME
    return
end

try
    [vdep, varr, ~, vp2, r1, ~] = leg_vel(a, b, t1, t2, dir, bodies, mu, ys);

    rp_AU = perihelion_AU(r1, vdep, mu, AU);
    if rp_AU < opts.MIN_PERIHELION_AU
        return
    end

    arr_vinf_last = varr(:) - vp2(:);
    node = pack_node([a b], [t1 t2], dir, arr_vinf_last, NaN, NaN, rp_AU);
    node = update_node_scores_and_rank(node, cfg, opts, bodies, mu, ys);

catch ME
    fail_msg = sprintf('SEED FAIL %d -> %d | TOF %.2f | dir %d | %s', ...
        a, b, tof, dir, ME.message);
end

end

%% =====================================================================
function children = expand_frontier(frontier, cfg, opts, bodies, mu, ys, AU)
%EXPAND_FRONTIER Expand every parent node by one additional flyby.
% Parallel version: each parent expansion is independent. Every worker builds
% a local child array; the client concatenates the cell array afterward.

if isempty(frontier)
    children = [];
    return
end

child_cells = cell(1, numel(frontier));

if opts.USE_PARALLEL && numel(frontier) > 1
    parfor ip = 1:numel(frontier)
        child_cells{ip} = expand_one_parent(frontier(ip), cfg, opts, bodies, mu, ys, AU);
    end
else
    for ip = 1:numel(frontier)
        child_cells{ip} = expand_one_parent(frontier(ip), cfg, opts, bodies, mu, ys, AU);
    end
end

children = concat_node_cells(child_cells);

end

%% =====================================================================
function children = expand_one_parent(parent, cfg, opts, bodies, mu, ys, AU)
%EXPAND_ONE_PARENT Worker-safe helper that expands one parent node.

children = [];
last_id = parent.seq(end);

for target = cfg.BODY_IDS
    if target == last_id
        continue
    end

    if count_body(parent.seq, target) >= opts.MAX_FLYBYS_PER_BODY
        continue
    end

    if target == 1000 && cfg.ALLOW_NYXAR_TERMINAL_ONLY && ~cfg.ALLOW_NYXAR_INTERIOR_PASS
        % In raw runs Nyxar is not useful as an interior flyby.
        % It can still appear only if the run ends there; beam cannot know that now.
        continue
    end

    for dir = cfg.DIRS
        cand = refine_next_tof_multi(parent, target, dir, cfg, opts, bodies, mu, ys, AU);

        for ic = 1:numel(cand)
            tof = cand(ic).tof;
            t1 = parent.t(end);
            t2 = t1 + tof;

            if t2 > cfg.MAX_TIME
                continue
            end

            if ~repeated_timing_ok(parent.seq, parent.t, target, t2, bodies, opts)
                continue
            end

            try
                [vdep, varr, vp1, vp2, r1, ~] = leg_vel(last_id, target, t1, t2, dir, bodies, mu, ys);

                rp_AU = perihelion_AU(r1, vdep, mu, AU);
                if rp_AU < opts.MIN_PERIHELION_AU
                    continue
                end

                vinf_in  = parent.arr_vinf_last(:);
                vinf_out = vdep(:) - vp1(:);

                [okfly, hR, delta_deg] = flyby_ok(last_id, vinf_in, vinf_out, bodies, opts, cfg);
                if ~okfly
                    continue
                end

                arr_vinf_last = varr(:) - vp2(:);
                child = pack_node([parent.seq target], [parent.t t2], [parent.dirs dir], ...
                    arr_vinf_last, hR, delta_deg, rp_AU);
                child = update_node_scores_and_rank(child, cfg, opts, bodies, mu, ys);

                children = append_node(children, child);

            catch
                % Keep the search clean. Bad Lambert arcs die here.
            end
        end
    end
end

end

%% =====================================================================
function cand = refine_next_tof_multi(parent, target, dir, cfg, opts, bodies, mu, ys, AU)
%REFINE_NEXT_TOF_MULTI Find one or more promising TOFs for the next leg.
% The method samples a coarse TOF grid, detects local minima of the flyby
% closure cost, then refines each selected interval with fminbnd(). This is
% more robust than running fminbnd() on one large interval.

cand = empty_candidate();

remaining = cfg.MAX_TIME - parent.t(end);
tof_hi = min(opts.TOF_MAX, remaining);
tof_lo = opts.TOF_MIN;

if tof_hi <= tof_lo
    return
end

cut = min(opts.TOF_SHORT_CUT_YR, tof_hi);
g1 = linspace(tof_lo, cut, opts.TOF_GRID_SHORT_N);
if tof_hi > cut + 1e-9
    g2 = linspace(cut, tof_hi, opts.TOF_GRID_LONG_N);
    grid = unique([g1 g2]);
else
    grid = unique(g1);
end

cost = inf(size(grid));
meta = repmat(empty_meta(), size(grid));

for k = 1:numel(grid)
    [cost(k), meta(k)] = closure_cost_at_tof(parent, target, dir, grid(k), cfg, opts, bodies, mu, ys, AU);
end

finite_idx = find(isfinite(cost));
if isempty(finite_idx)
    return
end

idxs = pick_best_grid_indices(cost, opts.N_TOP_INTERVALS);

out = empty_candidate();
for ii = 1:numel(idxs)
    idx = idxs(ii);
    lo = grid(max(1, idx-1));
    hi = grid(min(numel(grid), idx+1));

    if hi <= lo
        continue
    end

    try
        f = @(tof) closure_cost_scalar(parent, target, dir, tof, cfg, opts, bodies, mu, ys, AU);
        [tof_star, cost_star] = fminbnd(f, lo, hi, optimset('TolX', 1e-8, 'Display', 'off'));
    catch
        tof_star = grid(idx);
        cost_star = cost(idx);
    end

    [~, met] = closure_cost_at_tof(parent, target, dir, tof_star, cfg, opts, bodies, mu, ys, AU);

    if met.miss_mag <= opts.VINFTOL * opts.SOFT_VINFTOL_FACTOR
        c = make_candidate(tof_star, cost_star, met);
        out = append_candidate(out, c);
    end
end

if isempty(out)
    return
end

% remove near-duplicates by TOF
out = unique_candidates_by_tof(out, 1e-5);
[~, ord] = sort([out.cost], 'ascend');
out = out(ord);
cand = out(1:min(opts.N_TOP_TOF, numel(out)));

end

%% =====================================================================
function y = closure_cost_scalar(parent, target, dir, tof, cfg, opts, bodies, mu, ys, AU)
[y, ~] = closure_cost_at_tof(parent, target, dir, tof, cfg, opts, bodies, mu, ys, AU);
end

%% =====================================================================
function [cost, met] = closure_cost_at_tof(parent, target, dir, tof, cfg, opts, bodies, mu, ys, AU)
%CLOSURE_COST_AT_TOF Objective used during TOF refinement.
% Main term: absolute mismatch between outgoing and incoming V-infinity
% magnitudes at the flyby body. Soft terms penalize impossible turn angles and
% very long TOFs, guiding the optimizer toward useful physical closures.

met = empty_meta();
cost = Inf;

try
    t1 = parent.t(end);
    t2 = t1 + tof;

    if t2 > cfg.MAX_TIME || tof <= 0
        return
    end

    [vdep, ~, vp1, ~, r1, ~] = leg_vel(parent.seq(end), target, t1, t2, dir, bodies, mu, ys);

    rp_AU = perihelion_AU(r1, vdep, mu, AU);
    if rp_AU < opts.MIN_PERIHELION_AU
        return
    end

    vinf_in  = parent.arr_vinf_last(:);
    vinf_out = vdep(:) - vp1(:);

    nin = norm(vinf_in);
    nout = norm(vinf_out);
    if nin <= 0 || nout <= 0
        return
    end

    miss_mag = abs(nout - nin);
    delta_deg = turn_angle_deg(vinf_in, vinf_out);

    [turn_violation_deg, ~, ~] = turn_violation(parent.seq(end), vinf_in, vinf_out, bodies, opts, cfg);

    % Main term remains Vinf magnitude closure.
    % Small soft penalties guide fminbnd toward physically useful roots.
    cost = miss_mag ...
         + opts.COST_TURN_W * turn_violation_deg^2 ...
         + opts.COST_TOF_W  * max(0, tof - opts.TARGET_RECENT_TOF);

    met.miss_mag = miss_mag;
    met.delta_deg = delta_deg;
    met.turn_violation_deg = turn_violation_deg;
    met.rp_AU = rp_AU;

    if ~isfinite(cost)
        cost = Inf;
    end
catch
    cost = Inf;
end

end

%% =====================================================================
function idxs = pick_best_grid_indices(cost, N)

idxs = [];
finite_idx = find(isfinite(cost));
if isempty(finite_idx)
    return
end

% local minima first
for k = finite_idx(:).'
    left_ok  = (k == 1) || cost(k) <= cost(k-1);
    right_ok = (k == numel(cost)) || cost(k) <= cost(k+1);
    if left_ok && right_ok
        idxs(end+1) = k; %#ok<AGROW>
    end
end

if isempty(idxs)
    [~, ord] = sort(cost(finite_idx), 'ascend');
    idxs = finite_idx(ord(1:min(N, numel(ord))));
else
    [~, ord] = sort(cost(idxs), 'ascend');
    idxs = idxs(ord(1:min(N, numel(ord))));
end

end

%% =====================================================================
function [ok, h_over_R, delta_deg] = flyby_ok(id, vinf_in, vinf_out, bodies, opts, cfg)
%FLYBY_OK Strict patched-conic flyby validation.
% Checks V-infinity magnitude equality, required turn angle, and pericenter
% altitude. For a massless Nyxar flyby, vector continuity is required.

ok = false;
h_over_R = NaN;
delta_deg = NaN;

b = bodies.(sprintf('id_%d', id));
vin = norm(vinf_in);
vout = norm(vinf_out);

if vin <= 0 || vout <= 0
    return
end

if abs(vin - vout) > opts.VINFTOL
    return
end

if b.mu <= 0 || b.radius <= 0
    if cfg.ALLOW_NYXAR_INTERIOR_PASS && id == 1000
        dir_err = norm(vinf_out(:) - vinf_in(:));
        ok = dir_err <= opts.VINFTOL;
        h_over_R = Inf;
        delta_deg = 0;
    end
    return
end

c = dot(vinf_in, vinf_out)/(vin*vout);
c = max(-1, min(1, c));
delta = acos(c);
delta_deg = rad2deg(delta);

if delta < 1e-12
    return
end

rp = b.mu/vin^2 * (1/sin(delta/2) - 1);
h_over_R = rp/b.radius - 1;

ok = isfinite(h_over_R) && ...
     h_over_R >= opts.MIN_H_OVER_R && ...
     h_over_R <= opts.MAX_H_OVER_R;

end

%% =====================================================================
function [viol_deg, delta_min_deg, delta_max_deg] = turn_violation(id, vinf_in, vinf_out, bodies, opts, cfg)

viol_deg = 0;
delta_min_deg = 0;
delta_max_deg = 180;

b = bodies.(sprintf('id_%d', id));

if b.mu <= 0 || b.radius <= 0
    if cfg.ALLOW_NYXAR_INTERIOR_PASS && id == 1000
        viol_deg = turn_angle_deg(vinf_in, vinf_out);
        delta_min_deg = 0;
        delta_max_deg = 0;
    end
    return
end

vin = 0.5*(norm(vinf_in) + norm(vinf_out));
if vin <= 0
    viol_deg = Inf;
    return
end

rp_min = b.radius * (1 + opts.MIN_H_OVER_R);
rp_max = b.radius * (1 + opts.MAX_H_OVER_R);

delta_max = 2*asin(1/(1 + rp_min*vin^2/b.mu));
delta_min = 2*asin(1/(1 + rp_max*vin^2/b.mu));

delta_min_deg = rad2deg(delta_min);
delta_max_deg = rad2deg(delta_max);
delta_req_deg = turn_angle_deg(vinf_in, vinf_out);

if delta_req_deg < delta_min_deg
    viol_deg = delta_min_deg - delta_req_deg;
elseif delta_req_deg > delta_max_deg
    viol_deg = delta_req_deg - delta_max_deg;
else
    viol_deg = 0;
end

end

%% =====================================================================
function node = update_node_scores_and_rank(node, cfg, opts, bodies, mu, ys)

[J_raw, J_official, info] = score_node_dual(node, cfg, bodies, mu, ys);
node.J_raw = J_raw;
node.J_official = J_official;
node.score_info = info;

if cfg.USE_GRAND_TOUR_BONUS
    node.J = node.J_official;
else
    node.J = node.J_raw;
end

node.rank = rank_node_smart(node, cfg, opts, bodies);
node.score_rate = node.J_raw / max(node.t(end), 1e-9);

end

%% =====================================================================
function [J_raw, J_official, info] = score_node_dual(node, cfg, bodies, mu, ys)
%SCORE_NODE_DUAL Compute raw and optional official scores for one node.
% J_raw ignores the Grand Tour bonus. J_official applies the 1.2 multiplier
% only if the configuration enables it and all bodies are present.

seq = node.seq;
t = node.t;
dirs = node.dirs;

n = numel(seq);
nleg = n - 1;

R = zeros(n,3);
VP = zeros(n,3);

for k = 1:n
    s = KM_solver(t(k)*ys, bodies.(sprintf('id_%d', seq(k))));
    R(k,:)  = s(1:3).';
    VP(k,:) = s(4:6).';
end

vdep = zeros(nleg,3);
varr = zeros(nleg,3);

for leg = 1:nleg
    [vd, va] = leg_vel(seq(leg), seq(leg+1), t(leg), t(leg+1), dirs(leg), bodies, mu, ys);
    vdep(leg,:) = vd(:).';
    varr(leg,:) = va(:).';
end

vinf = zeros(1,n);
for k = 1:n
    if k == 1
        vinf(k) = norm(vdep(1,:) - VP(k,:));
    elseif k <= nleg
        vin = norm(varr(k-1,:) - VP(k,:));
        vou = norm(vdep(k,:)   - VP(k,:));
        vinf(k) = 0.5*(vin + vou);
    else
        vinf(k) = norm(varr(k-1,:) - VP(k,:));
    end
end

rhat = R ./ vecnorm(R,2,2);
all_ids = [4 5 6 7 8 9 10 1000];

J_raw = 0;
by_body = struct();

for id = all_ids
    idx = find(seq == id);
    if isempty(idx)
        continue
    end

    idx = idx(1:min(13,numel(idx)));
    b = bodies.(sprintf('id_%d', id));
    w = b.weight;

    prev_dirs = [];
    contrib = 0;

    for ii = 1:numel(idx)
        k = idx(ii);
        V = vinf(k);

        F = 0.2 + exp(-V/13)/(1 + exp(-5*(V - 1.5)));

        if ii == 1
            S = 1.0;
        else
            dots = prev_dirs * rhat(k,:)';
            dots = max(-1,min(1,dots));
            ang = acosd(dots);
            S = 0.1 + 0.9/(1 + 10*sum(exp(-(ang.^2)/50)));
        end

        dJ = w*S*F;
        contrib = contrib + dJ;
        J_raw = J_raw + dJ;
        prev_dirs = [prev_dirs; rhat(k,:)]; %#ok<AGROW>
    end

    by_body.(sprintf('id_%d', id)) = contrib;
end

visited = unique(seq);
grand_tour_hit = all(ismember(all_ids, visited));
bonus = 1.0;
if cfg.USE_GRAND_TOUR_BONUS && grand_tour_hit
    bonus = 1.2;
end
J_official = J_raw * bonus;

info.vinf = vinf;
info.grand_tour_hit = grand_tour_hit;
info.bonus = bonus;
info.by_body = by_body;

end

%% =====================================================================
function R = rank_node_smart(node, cfg, opts, bodies)
%RANK_NODE_SMART Heuristic score used for beam pruning.
%
% IMPORTANT:
% The final objective is still J_raw. This rank is only a heuristic used to
% decide which partial trajectories survive. The weights are defined in
% default_opts(), section F: "RANKING PARAMETERS / HEURISTIC WEIGHTS".

J = node.J_raw;
t_now = node.t(end);
n_vis = numel(node.seq);

yrs_left = max(0, cfg.MAX_TIME - t_now);
slots_left = max(0, cfg.MAX_DEPTH - n_vis);

if yrs_left <= 0 || slots_left <= 0
    R = opts.RANK_WEIGHT_RAW_SCORE * J - 1000*max(0,t_now-cfg.MAX_TIME);
    return
end

score_rate = J / max(t_now, 1e-6);

if n_vis >= 2
    dt = diff(node.t);
    avg_tof = mean(dt);
    recent_tof = mean(dt(max(1,end-2):end));
else
    avg_tof = opts.TARGET_MEAN_TOF;
    recent_tof = opts.TARGET_MEAN_TOF;
end

% A path at high mission time with too few flybys is penalized.
target_vis_now = 1 + floor(t_now / opts.TARGET_MEAN_TOF);
lag = max(0, target_vis_now - n_vis);

% Estimate realistic remaining useful slots.
slots_by_time = max(0, floor(yrs_left / max(5, opts.TARGET_RECENT_TOF)) + 1);
slots_realistic = min(slots_left, slots_by_time);

future_gains = [];
for id = cfg.BODY_IDS
    n_used = sum(node.seq == id);
    n_left_rule = max(0, opts.MAX_FLYBYS_PER_BODY - n_used);
    if n_left_rule <= 0
        continue
    end

    b = bodies.(sprintf('id_%d', id));
    w = b.weight;
    T = b.period_year;

    gap = max(opts.MIN_INSERT_GAP, opts.REPEAT_FRAC*T);
    if n_used == 0
        n_time_real = 1 + floor(max(0, yrs_left - opts.MIN_INSERT_GAP)/gap);
        season0 = 1.00;
    elseif n_used == 1
        n_time_real = floor(yrs_left/gap);
        season0 = 0.90;
    elseif n_used == 2
        n_time_real = floor(yrs_left/gap);
        season0 = 0.72;
    else
        n_time_real = floor(yrs_left/gap);
        season0 = max(0.35, 0.72*0.85^(n_used-2));
    end

    n_real = min([n_left_rule, n_time_real, slots_realistic]);
    if n_real <= 0
        continue
    end

    base = w * opts.RANK_F_REF;
    if ismember(id, opts.HEAVY_IDS)
        base = 1.10*base;
    end
    if ismember(id, opts.LOWVALUE_IDS)
        base = 0.70*base;
    end

    for k = 1:n_real
        future_gains(end+1) = base * max(0.25, season0*0.88^(k-1)); %#ok<AGROW>
    end
end

future_gains = sort(future_gains,'descend');
future_bonus = opts.RANK_FUTURE_SCALE * sum(future_gains(1:min(slots_realistic,numel(future_gains))));

heavy_unique = numel(intersect(unique(node.seq), opts.HEAVY_IDS));
all_unique = numel(unique(node.seq));

slow_pen   = opts.PEN_SLOW_W * max(0, avg_tof    - opts.TARGET_MEAN_TOF)^2;
recent_pen = opts.PEN_RECENT_SLOW_W * max(0, recent_tof - opts.TARGET_RECENT_TOF)^2;

pen_low_4 = opts.PEN_LOW_4_W * max(0, sum(node.seq == 4) - 3);
pen_low_5 = 0.4 * opts.PEN_LOW_5_W * max(0, sum(node.seq == 5) - 8);
pen_nyxar = opts.PEN_NYXAR_W * sum(node.seq == 1000) * (~cfg.USE_GRAND_TOUR_BONUS);

seq = node.seq;
n9 = sum(seq == 9);
n8 = sum(seq == 8);
n4 = sum(seq == 4);

bonus_family249 = 0;

if ~isempty(seq) && seq(1) == 9
    bonus_family249 = bonus_family249 + 70;
end

if numel(seq) >= 2 && seq(1) == 9 && seq(2) == 8
    bonus_family249 = bonus_family249 + 60;
end

bonus_family249 = bonus_family249 + 30 * min(n8, 3);

if n9 >= 2
    bonus_family249 = bonus_family249 + 45;
end

if ~isempty(seq) && ~isempty(node.t) && seq(end) == 9 && node.t(end) > 260
    bonus_family249 = bonus_family249 + 70;
end

if numel(seq) >= 2 && seq(end-1) == 8 && seq(end) == 9
    bonus_family249 = bonus_family249 + 60;
end

pen_family249 = 0;

if n9 == 0
    pen_family249 = pen_family249 + 120;
end

if n8 == 0
    pen_family249 = pen_family249 + 80;
end

if ~isempty(node.t) && node.t(end) > 260 && n9 < 2
    pen_family249 = pen_family249 + 80;
end

if n4 > 3
    pen_family249 = pen_family249 + 25*(n4-3);
end

late_pen = opts.PEN_LATE_W * (t_now / cfg.MAX_TIME)^4;
lag_pen  = opts.PEN_LAG_W * lag^2;

R = opts.RANK_WEIGHT_RAW_SCORE * J ...
  + opts.RANK_WEIGHT_SCORE_RATE * score_rate ...
  + future_bonus ...
  + opts.RANK_WEIGHT_HEAVY_UNIQUE * heavy_unique ...
  + opts.RANK_WEIGHT_ALL_UNIQUE * all_unique ...
  - lag_pen ...
  - slow_pen ...
  - recent_pen ...
  - pen_low_4 ...
  - pen_low_5 ...
  - pen_nyxar ...
  - late_pen;

R = R + bonus_family249 - pen_family249;

end


%% =====================================================================
function out = prune_frontier_diverse(nodes, N, cfg, opts)
%PRUNE_FRONTIER_DIVERSE Keep the beam width under control without collapsing
% to one repeated family. Nodes are sorted by rank, but bucket protection keeps
% representatives with different last bodies, time bins, V-infinity bins, and
% heavy-body coverage.

if isempty(nodes)
    out = nodes;
    return
end

nodes = sort_nodes(nodes, 'rank');
N = min(N, numel(nodes));

% ---------------------------------------------------------------------
% RAW SCORE PROTECTION
% ---------------------------------------------------------------------
% This is the key score-hunting modification. A fraction of the beam is
% protected according to J_raw before bucket/rank pruning is applied.
% Tune opts.RAW_PROTECT_FRAC and opts.RAW_PROTECT_MIN in default_opts().
raw_protect_frac = opts.RAW_PROTECT_FRAC;
n_raw_keep = min(numel(nodes), max(opts.RAW_PROTECT_MIN, round(raw_protect_frac * N)));
[~, raw_ord] = sort([nodes.J_raw], 'descend');
raw_keep_idx = raw_ord(1:n_raw_keep);

% Bucket protection: last body + time + Vinf + heavy coverage.
keys = cell(1,numel(nodes));
for i = 1:numel(nodes)
    last = nodes(i).seq(end);
    tbin = floor(nodes(i).t(end) / opts.TIME_BIN_YR);
    vbin = floor(norm(nodes(i).arr_vinf_last) / opts.VINF_BIN);
    hcnt = numel(intersect(unique(nodes(i).seq), opts.HEAVY_IDS));
    keys{i} = sprintf('%d_%d_%d_%d', last, tbin, vbin, hcnt);
end

[uk,~,ic] = unique(keys);

% Start from the raw-score-protected nodes, then add bucket representatives.
keep_idx = raw_keep_idx(:).';

for k = 1:numel(uk)
    idx = find(ic == k);
    idx = idx(:).';
    [~, loc_ord] = sort([nodes(idx).rank], 'descend');
    idx = idx(loc_ord);
    nk = min(opts.KEEP_PER_BUCKET, numel(idx));
    keep_idx = [keep_idx idx(1:nk)]; %#ok<AGROW>
end

keep_idx = unique(keep_idx(:).', 'stable');

if numel(keep_idx) > N
    [~, ord] = sort([nodes(keep_idx).rank], 'descend');
    keep_idx = keep_idx(ord(1:N));
    out = nodes(keep_idx);
    return
end

all_idx = 1:numel(nodes);
rest_idx = setdiff(all_idx, keep_idx, 'stable');
rest_idx = rest_idx(:).';

if isfield(cfg, 'ELITE_FRAC')
    elite_frac = cfg.ELITE_FRAC;
else
    elite_frac = opts.ELITE_FRAC;
end

if isfield(cfg, 'SOFTMAX_T')
    softmax_T = cfg.SOFTMAX_T;
else
    softmax_T = opts.SOFTMAX_T;
end

if cfg.USE_STOCHASTIC_BEAM && ~isempty(rest_idx)
    n_left = N - numel(keep_idx);
    n_elite = min(n_left, max(0, round(elite_frac * N) - numel(keep_idx)));

    elite_rest = rest_idx(1:min(n_elite, numel(rest_idx)));
        elite_rest = elite_rest(:).';
        pool = setdiff(rest_idx, elite_rest, 'stable');
        pool = pool(:).';

    n_sample = N - numel(keep_idx) - numel(elite_rest);
    if n_sample > 0 && ~isempty(pool)
        ranks = [nodes(pool).rank];
        if ~isfinite(softmax_T) || softmax_T <= 0
            softmax_T = 1.0;
        end
        ranks(~isfinite(ranks)) = -Inf;

        n_pick = min(n_sample, numel(pool));
        if all(~isfinite(ranks))
            pick_local = randperm(numel(pool), n_pick);
        else
            rmax = max(ranks(isfinite(ranks)));
            z = (ranks - rmax) / softmax_T;
            z(~isfinite(z)) = -Inf;
            p = exp(z);
            p(~isfinite(p) | p < 0) = 0;
            ps = sum(p);
            if ps <= 0
                pick_local = randperm(numel(pool), n_pick);
            else
                p = p / ps;
                pick_local = weighted_sample_no_replace(p, n_pick);
            end
        end
        sampled = pool(pick_local);
        sampled = sampled(:).';
    else
        sampled = zeros(1,0);
    end

    final_idx = [keep_idx(:).' elite_rest(:).' sampled(:).'];
else
    final_idx = [keep_idx(:).' rest_idx(:).'];
end

final_idx = unique(final_idx(:).', 'stable');
final_idx = final_idx(1:min(N,numel(final_idx)));
out = nodes(final_idx);

end

%% =====================================================================
function pick = weighted_sample_no_replace(p, n)

p = p(:).';
n = min(max(0, n), numel(p));
if n == 0
    pick = zeros(1,0);
    return
end

pick = zeros(1,n);
avail = 1:numel(p);
w = p;
w(~isfinite(w) | w < 0) = 0;
if sum(w) <= 0
    w = ones(size(w));
end

for k = 1:n
    sw = sum(w);
    if sw <= 0 || ~isfinite(sw)
        j = randi(numel(avail));
    else
        w = w / sw;
        r = rand;
        c = cumsum(w);
        j = find(r <= c, 1, 'first');
        if isempty(j)
            j = numel(avail);
        end
    end
    pick(k) = avail(j);
    avail(j) = [];
    w(j) = [];
    if isempty(avail)
        pick = pick(1:k);
        return
    end
end

end

%% =====================================================================
function ref_best = refine_sequence_times(node, cfg, opts, bodies, mu, ys, AU)
%REFINE_SEQUENCE_TIMES Local continuous refinement of a fixed sequence.
% This does not change the body order or Lambert directions. It perturbs the
% leg durations and keeps changes only if the fully validated node improves.

ref_best = node;
dt0 = diff(node.t);

if isempty(dt0) || numel(dt0) ~= numel(node.dirs)
    return
end

for rr = 1:opts.REFINE_RESTARTS
    if rr == 1
        dt = dt0;
    else
        dt = dt0 .* (1 + opts.REFINE_STEP_FRAC*randn(size(dt0)));
        dt = max(opts.TOF_MIN, dt);
        if sum(dt) > cfg.MAX_TIME
            dt = dt * (0.995*cfg.MAX_TIME/sum(dt));
        end
    end

    cand = try_build_node_from_times(node.seq, [0 cumsum(dt)], node.dirs, cfg, opts, bodies, mu, ys, AU);
    if isempty(cand)
        continue
    end

    step = max(opts.REFINE_MIN_STEP_YR, opts.REFINE_STEP_FRAC*dt);

    for sweep = 1:opts.REFINE_SWEEPS
        improved = false;
        for j = 1:numel(dt)
            trials = [dt(j)-step(j), dt(j)+step(j)];
            for tv = trials
                if tv < opts.TOF_MIN
                    continue
                end
                dt_try = dt;
                dt_try(j) = tv;
                if sum(dt_try) > cfg.MAX_TIME
                    continue
                end
                cand_try = try_build_node_from_times(node.seq, [0 cumsum(dt_try)], node.dirs, cfg, opts, bodies, mu, ys, AU);
                if ~isempty(cand_try) && cand_try.J_raw > cand.J_raw
                    cand = cand_try;
                    dt = dt_try;
                    improved = true;
                end
            end
        end
        step = 0.5*step;
        if ~improved
            % keep shrinking step
        end
    end

    if cand.J_raw > ref_best.J_raw
        ref_best = cand;
    end
end

end

%% =====================================================================
function node = try_build_node_from_times(seq, t, dirs, cfg, opts, bodies, mu, ys, AU)
%TRY_BUILD_NODE_FROM_TIMES Rebuild and validate a full node from fixed times.
% Used by the local refinement step to ensure that every flyby remains
% dynamically consistent after timing perturbations.

node = [];

if any(diff(t) <= 0) || t(end) > cfg.MAX_TIME
    return
end

n = numel(seq);
if numel(dirs) ~= n-1
    return
end

% Check repeated timing.
for k = 2:n
    if ~repeated_timing_ok(seq(1:k-1), t(1:k-1), seq(k), t(k), bodies, opts)
        return
    end
end

arr_vinf_last = [];
hR_last = NaN;
delta_last = NaN;
rp_AU_last = NaN;

% First leg.
try
    [~, varr, ~, vp2, r1, ~] = leg_vel(seq(1), seq(2), t(1), t(2), dirs(1), bodies, mu, ys);
    rp_AU = perihelion_AU(r1, varr, mu, AU); %#ok<NASGU>
catch
    return
end

% Full leg/flyby validation.
for leg = 1:n-1
    try
        [vdep, varr, vp1, vp2, r1, ~] = leg_vel(seq(leg), seq(leg+1), t(leg), t(leg+1), dirs(leg), bodies, mu, ys);
    catch
        return
    end

    rp_AU_leg = perihelion_AU(r1, vdep, mu, AU);
    if rp_AU_leg < opts.MIN_PERIHELION_AU
        return
    end

    if leg >= 2
        % Validate flyby at seq(leg).
        [~, varr_prev, ~, vp_here, ~, ~] = leg_vel(seq(leg-1), seq(leg), t(leg-1), t(leg), dirs(leg-1), bodies, mu, ys);
        vinf_in = varr_prev(:) - vp_here(:);
        vinf_out = vdep(:) - vp1(:);
        [okfly, hR_last, delta_last] = flyby_ok(seq(leg), vinf_in, vinf_out, bodies, opts, cfg);
        if ~okfly
            return
        end
    end

    arr_vinf_last = varr(:) - vp2(:);
    rp_AU_last = rp_AU_leg;
end

node = pack_node(seq, t, dirs, arr_vinf_last, hR_last, delta_last, rp_AU_last);
node = update_node_scores_and_rank(node, cfg, opts, bodies, mu, ys);

end

%% =====================================================================
function [vdep, varr, vp1, vp2, r1, r2] = leg_vel(id1, id2, t1, t2, dir, bodies, mu, ys)
%LEG_VEL Convenience wrapper around ephemerides + Lambert solver.
% Returns departure/arrival spacecraft velocities and body velocities at both
% event times, all as column vectors.

s1 = KM_solver(t1*ys, bodies.(sprintf('id_%d', id1)));
s2 = KM_solver(t2*ys, bodies.(sprintf('id_%d', id2)));

r1 = s1(1:3);
r2 = s2(1:3);
vp1 = s1(4:6);
vp2 = s2(4:6);

[vdep, varr] = lambert_universal(mu, t1*ys, r1, t2*ys, r2, dir);

vdep = real(vdep(:));
varr = real(varr(:));
vp1  = vp1(:);
vp2  = vp2(:);
r1   = r1(:);
r2   = r2(:);

if any(~isfinite(vdep)) || any(~isfinite(varr))
    error('Lambert returned non-finite velocity')
end

end

%% =====================================================================
function rp_AU = perihelion_AU(r, v, mu, AU)
%PERIHELION_AU Estimate the heliocentric perihelion radius of a Lambert arc.
% Used as a fast screening constraint before accepting a candidate leg.

r = r(:);
v = v(:);

h = cross(r,v);
h2 = dot(h,h);

evec = cross(v,h)/mu - r/norm(r);
e = norm(evec);

p = h2/mu;

if e < 1e-12
    rp = p;
else
    rp = p/(1+e);
end

rp_AU = rp/AU;

if ~isfinite(rp_AU) || rp_AU <= 0
    rp_AU = Inf;
end

end

%% =====================================================================
function ok = repeated_timing_ok(seq, t, target, t_new, bodies, opts)
%REPEATED_TIMING_OK Enforce the minimum interval for repeated flybys.
% If the same target was visited before, the new visit must occur at least
% one third of that body's orbital period later.

ok = true;
idx = find(seq == target, 1, 'last');
if isempty(idx)
    return
end

Tbody = bodies.(sprintf('id_%d',target)).period_year;
min_gap = opts.REPEAT_FRAC * Tbody;

if t_new - t(idx) < min_gap
    ok = false;
end

end

%% =====================================================================
function n = count_body(seq, id)
n = sum(seq == id);
end

%% =====================================================================
function node = pack_node(seq, t, dirs, arr_vinf_last, hR_last, delta_last, rp_AU_last)

node = struct();
node.seq = seq;
node.t = t;
node.dirs = dirs;
node.arr_vinf_last = arr_vinf_last(:);
node.J = NaN;
node.J_raw = NaN;
node.J_official = NaN;
node.rank = NaN;
node.score_rate = NaN;
node.hR_last = hR_last;
node.delta_last = delta_last;
node.rp_AU_last = rp_AU_last;
node.ok = true;
node.score_info = struct();

end

%% =====================================================================
function out = append_node(out, node)
if isempty(out)
    out = node;
else
    out(end+1) = node; %#ok<AGROW>
end
end


%% =====================================================================
function out = concat_node_cells(cells_in)
%CONCAT_NODE_CELLS Concatenate a cell array of node arrays.
% This helper keeps parfor loops valid because workers write to sliced cells
% instead of growing one shared struct array.

out = [];
if isempty(cells_in)
    return
end

for i = 1:numel(cells_in)
    item = cells_in{i};
    if isempty(item)
        continue
    end
    if isempty(out)
        out = item;
    else
        out = [out item]; %#ok<AGROW>
    end
end

end


%% =====================================================================
function c = make_candidate(tof, cost, met)
c = struct();
c.tof = tof;
c.cost = cost;
c.miss_mag = met.miss_mag;
c.delta_deg = met.delta_deg;
c.turn_violation_deg = met.turn_violation_deg;
c.rp_AU = met.rp_AU;
end

%% =====================================================================
function c = empty_candidate()
c = struct('tof', {}, 'cost', {}, 'miss_mag', {}, 'delta_deg', {}, ...
           'turn_violation_deg', {}, 'rp_AU', {});
end

%% =====================================================================
function c = append_candidate(c, item)
if isempty(c)
    c = item;
else
    c(end+1) = item; %#ok<AGROW>
end
end

%% =====================================================================
function m = empty_meta()
m = struct();
m.miss_mag = Inf;
m.delta_deg = NaN;
m.turn_violation_deg = Inf;
m.rp_AU = Inf;
end

%% =====================================================================
function out = unique_candidates_by_tof(cand, tol)

out = empty_candidate();
for i = 1:numel(cand)
    duplicate = false;
    for j = 1:numel(out)
        if abs(cand(i).tof - out(j).tof) < tol
            duplicate = true;
            if cand(i).cost < out(j).cost
                out(j) = cand(i);
            end
            break
        end
    end
    if ~duplicate
        out = append_candidate(out, cand(i));
    end
end

end

%% =====================================================================
function nodes = sort_nodes(nodes, mode)

if isempty(nodes)
    return
end

switch lower(mode)
    case 'rank'
        q = [nodes.rank];
    case 'raw'
        q = [nodes.J_raw];
    case 'official'
        q = [nodes.J_official];
    otherwise
        q = [nodes.J];
end

[~, idx] = sort(q, 'descend');
nodes = nodes(idx);

end

%% =====================================================================
function archive = merge_archive(archive, new_nodes, maxN)

if isempty(new_nodes)
    return
end

if isempty(archive)
    archive = new_nodes;
else
    archive = [archive new_nodes]; %#ok<AGROW>
end

archive = prune_archive_unique(archive, maxN);

end

%% =====================================================================
function archive = prune_archive_unique(archive, maxN)

if isempty(archive)
    return
end

archive = sort_nodes(archive, 'raw');
seen = containers.Map('KeyType','char','ValueType','logical');
keep = false(1,numel(archive));

for i = 1:numel(archive)
    key = node_signature(archive(i));
    if ~isKey(seen, key)
        seen(key) = true;
        keep(i) = true;
    end
end

archive = archive(keep);
archive = archive(1:min(maxN,numel(archive)));

end

%% =====================================================================
function key = node_signature(node)
% Include sequence + rounded times + dirs. This avoids duplicate clones but
% still keeps useful timing variants.
seq_s = sprintf('%d-', node.seq);
t_s = sprintf('%.3f-', node.t);
d_s = sprintf('%d-', node.dirs);
key = [seq_s '|' t_s '|' d_s];
end

%% =====================================================================
function elites = take_top_unique(archive, N)

archive = prune_archive_unique(archive, max(5*N, N));
archive = sort_nodes(archive, 'raw');
elites = archive(1:min(N,numel(archive)));

end

%% =====================================================================
function save_run_result(res, cfg, opts)

if isempty(res.best)
    return
end

safe = matlab.lang.makeValidName(cfg.NAME);
matfile = fullfile(opts.OUTDIR, sprintf('run_%s.mat', safe));
save(matfile, 'res', 'cfg')

fprintf('Saved run result: %s\n', matfile)

end

%% =====================================================================
function write_archive_summary(archive, filename, bodies)

fid = fopen(filename, 'w');
if fid < 0
    warning('Could not write archive summary: %s', filename)
    return
end

fprintf(fid, 'rank,J_raw,J_official,t_end,nflybys,score_rate,sequence,names,times,dirs\n');

for i = 1:numel(archive)
    node = archive(i);
    names = cell(1,numel(node.seq));
    for k = 1:numel(node.seq)
        names{k} = body_label(node.seq(k), bodies);
    end

    seq_s = strtrim(sprintf('%d ', node.seq));
    names_s = strjoin(names, '->');
    t_s = strtrim(sprintf('%.6f ', node.t));
    d_s = strtrim(sprintf('%d ', node.dirs));

    fprintf(fid, '%d,%.12f,%.12f,%.9f,%d,%.12f,"%s","%s","%s","%s"\n', ...
        i, node.J_raw, node.J_official, node.t(end), numel(node.seq), ...
        node.J_raw/max(node.t(end),1e-9), seq_s, names_s, t_s, d_s);
end

fclose(fid);

end

%% =====================================================================
function export_top_csvs(archive, opts, N)
%EXPORT_TOP_CSVS Convert the best archived nodes to web-submission CSV files.
% The actual file formatting and any calibrated position patching are handled
% by solution_to_csv_patched().

if isempty(archive)
    return
end

archive = sort_nodes(archive, 'raw');
N = min(N, numel(archive));

csv_dir = fullfile(opts.OUTDIR, 'csv');
if ~isfolder(csv_dir)
    mkdir(csv_dir)
end

fprintf('\n=== EXPORT TOP %d CSV ===\n', N)

for i = 1:N
    node = archive(i);
    sol.sequence = node.seq;
    sol.t_flyby = node.t;
    sol.Lambert_dir = node.dirs;

    seq_tag = strjoin(arrayfun(@num2str, sol.sequence, 'UniformOutput', false), '-');
    outfile = fullfile(csv_dir, sprintf('rank%02d_J%.6f_%s.csv', i, node.J_raw, seq_tag));

    try
        solution_to_csv_patched(sol, outfile);
        fprintf('#%02d wrote %s\n', i, outfile)
    catch ME
        fprintf('#%02d export failed: %s\n', i, ME.message)
    end
end

end

%% =====================================================================
function print_best(nodes, N)

if isempty(nodes)
    fprintf('No nodes.\n')
    return
end

N = min(N, numel(nodes));
for k = 1:N
    fprintf('#%02d raw=%.6f off=%.6f rank=%.6f t=%.2f rate=%.4f seq=[%s]\n', ...
        k, nodes(k).J_raw, nodes(k).J_official, nodes(k).rank, nodes(k).t(end), ...
        nodes(k).J_raw/max(nodes(k).t(end),1e-9), sprintf('%d ', nodes(k).seq))
end

end

%% =====================================================================
function print_family249_diag(frontier)

if isempty(frontier)
    fprintf('Family249 diag | 2x9=0 | 3x8=0 | start98=0 | end89=0\n')
    return
end

frontier_has_2x9 = sum(arrayfun(@(n) sum(n.seq==9)>=2, frontier));
frontier_has_3x8 = sum(arrayfun(@(n) sum(n.seq==8)>=3, frontier));
frontier_start_98 = sum(arrayfun(@(n) numel(n.seq)>=2 && n.seq(1)==9 && n.seq(2)==8, frontier));
frontier_end_89 = sum(arrayfun(@(n) numel(n.seq)>=2 && n.seq(end-1)==8 && n.seq(end)==9, frontier));

fprintf('Family249 diag | 2x9=%d | 3x8=%d | start98=%d | end89=%d\n', ...
    frontier_has_2x9, frontier_has_3x8, frontier_start_98, frontier_end_89);

end

%% =====================================================================
function disp_best(node)

if isempty(node)
    fprintf('No best node.\n')
    return
end

fprintf('J_raw      = %.12f\n', node.J_raw)
fprintf('J_official = %.12f\n', node.J_official)
fprintf('rank       = %.12f\n', node.rank)
fprintf('t_end      = %.9f yr\n', node.t(end))
fprintf('score/year = %.9f\n', node.J_raw/max(node.t(end),1e-9))
fprintf('seq        = [%s]\n', sprintf('%d ', node.seq))
fprintf('times      = [%s]\n', sprintf('%.6f ', node.t))
fprintf('dirs       = [%s]\n', sprintf('%d ', node.dirs))

end

%% =====================================================================
function print_sequence_details(node, bodies, mu, ys)

if isempty(node)
    return
end

seq = node.seq;
t = node.t;
dirs = node.dirs;

n = numel(seq);
nleg = n - 1;

VP = zeros(n,3);
for k = 1:n
    s = KM_solver(t(k)*ys, bodies.(sprintf('id_%d',seq(k))));
    VP(k,:) = s(4:6).';
end

vdep = zeros(nleg,3);
varr = zeros(nleg,3);
for leg = 1:nleg
    [vd,va] = leg_vel(seq(leg),seq(leg+1),t(leg),t(leg+1),dirs(leg),bodies,mu,ys);
    vdep(leg,:) = vd(:).';
    varr(leg,:) = va(:).';
end

names = cell(1,n);
for k = 1:n
    names{k} = body_label(seq(k), bodies);
end

fprintf('\nSequence names:\n')
fprintf('%s\n', strjoin(names, ' -> '))

fprintf('\nDetails:\n')
for k = 1:n
    if k == 1
        vinf_in = norm(vdep(1,:) - VP(k,:));
        vinf_out = NaN;
        delta = NaN;
    elseif k <= nleg
        vinf_in_vec = varr(k-1,:) - VP(k,:);
        vinf_out_vec = vdep(k,:) - VP(k,:);
        vinf_in = norm(vinf_in_vec);
        vinf_out = norm(vinf_out_vec);
        delta = turn_angle_deg(vinf_in_vec, vinf_out_vec);
    else
        vinf_in = norm(varr(k-1,:) - VP(k,:));
        vinf_out = NaN;
        delta = NaN;
    end

    if k <= nleg
        tof = t(k+1) - t(k);
    else
        tof = NaN;
    end

    fprintf('%2d | %-12s | t=%9.4f | TOF=%9.4f | vinf_in=%9.4f | vinf_out=%9.4f | turn=%8.3f deg\n', ...
        k, names{k}, t(k), tof, vinf_in, vinf_out, delta)
end

end

%% =====================================================================
function delta = turn_angle_deg(vinf_in, vinf_out)

vin = norm(vinf_in);
vout = norm(vinf_out);

if vin <= 0 || vout <= 0
    delta = NaN;
    return
end

c = dot(vinf_in, vinf_out)/(vin*vout);
c = max(-1, min(1, c));
delta = rad2deg(acos(c));

end

%% =====================================================================
function name = body_label(id, bodies)

key = sprintf('id_%d', id);
if isfield(bodies, key) && isfield(bodies.(key), 'name')
    name = bodies.(key).name;
else
    name = key;
end

end

%% =====================================================================
function [root, paths] = setup_paths(script_dir)
%SETUP_PATHS Locate project folders and add them to the path.
% addpath(..., '-begin') is used so the intended project functions take
% priority over older copies elsewhere on the MATLAB path.

root = find_project_root_relaxed(script_dir);

paths = struct();
paths.script_dir = script_dir;
paths.flyby_dir = fullfile(root, 'flyby_functions');
paths.submission_dir = fullfile(root, 'submission');
paths.submission_utils_dir = fullfile(paths.submission_dir, 'utils');
paths.orbit_calibration_dir = fullfile(root, 'orbit_calibration');
paths.functions_dir = fullfile(root, 'functions');

add_if_dir(paths.flyby_dir)
add_if_dir(paths.submission_dir)
add_if_dir(paths.submission_utils_dir)
add_if_dir(paths.orbit_calibration_dir)
add_if_dir(paths.functions_dir)
add_if_dir(root)
add_if_dir(script_dir)

rehash toolboxcache

end

%% =====================================================================
function add_if_dir(p)
if isfolder(p)
    addpath(p, '-begin')
end
end

%% =====================================================================
function root = find_project_root_relaxed(start_dir)

root = start_dir;
for k = 1:10
    has_flyby = isfolder(fullfile(root,'flyby_functions')) || isfolder(fullfile(root,'functions','flyby_functions'));
    has_submission = isfolder(fullfile(root,'submission')) || isfile(fullfile(root,'solution_to_csv_patched.m'));
    has_bodies = isfile(fullfile(root,'bodies.mat')) || isfile(fullfile(root,'submission','bodies.mat'));

    if has_bodies && (has_flyby || has_submission)
        return
    end

    parent = fileparts(root);
    if strcmp(parent, root)
        break
    end
    root = parent;
end

% Fallback: current script folder.
root = start_dir;

end

%% =====================================================================
function bodies_file = find_bodies_file(root, paths, script_dir)

candidates = {
    fullfile(paths.submission_dir, 'bodies.mat')
    fullfile(paths.flyby_dir, 'bodies.mat')
    fullfile(root, 'bodies.mat')
    fullfile(script_dir, 'bodies.mat')
    };

for i = 1:numel(candidates)
    if isfile(candidates{i})
        bodies_file = candidates{i};
        return
    end
end

error('Missing bodies.mat. Put it near this script or in submission/flyby_functions.')

end
