%% ================== GENETIC ALGORITHM ==================

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
load(fullfile(thisDir, 'bodies.mat'), 'bodies');

y2sec = bodies.year_sec;

%% --- Problem setup ---

P1 = 8;   % flyby body ID
P2 = 9;   % next body ID

id_P1 = sprintf('id_%d', P1);
id_P2 = sprintf('id_%d', P2);

% Fixed initial conditions
t1 = 25;     % flyby arrival epoch [years]

state_P1 = KM_solver(t1 * y2sec, bodies.(id_P1));
state_P1 = state_P1(:);

v_planet = state_P1(4:6);
vinf_in  = [3.5; -1.5; 2.3];
v_in     = v_planet + vinf_in;

fprintf('Flyby at %s -> %s   |    t1 = %.0f years   |   v_inf = %.2f km/s\n', ...
    bodies.(id_P1).name, bodies.(id_P2).name, t1, norm(vinf_in));

%% --- Genetic Algorithm parameters ---

lb = [1, 0];        % [TOF_min, path_min]
ub = [80, 1];       % [TOF_max, path_max]

options = optimoptions('ga', ...
    'PopulationSize', 200, ...
    'MaxGenerations', 300, ...
    'FunctionTolerance', 1e-10, ...
    'EliteCount', 10, ...
    'Display', 'iter');   % this gives the 4-column GA output

fitness = @(x) flyby_cost_exercise(x, v_in, t1, id_P1, id_P2, bodies);

fprintf('\nRunning GA...\n');
tic;
[x_opt, fval] = ga(fitness, 2, [], [], [], [], lb, ub, [], [2], options);
elapsed = toc;

x_opt(2) = round(x_opt(2));
TOF  = x_opt(1);
path = x_opt(2);

fprintf('\nGA finished in %.1f seconds.\n', elapsed);
fprintf('TOF = %.4f years\n', TOF);
fprintf('Path = %d\n', path);
fprintf('Cost = %.4f\n', fval);

% Final physical check
[~, met] = flyby_cost_exercise(x_opt, v_in, t1, id_P1, id_P2, bodies);

fprintf('\n--- Flyby diagnostics ---\n');
fprintf('|vinf_in|      = %.9f km/s\n', met.vinf_in_mag);
fprintf('|vinf_out|     = %.9f km/s\n', met.vinf_out_mag);
fprintf('vinf mag error = %.9e km/s\n', met.mag_err);
fprintf('turn angle     = %.6f deg\n', met.delta_deg);
fprintf('h/R_body       = %.9f\n', met.h_over_R);
fprintf('|vinf_arr P2|  = %.9f km/s\n', met.vinf_arr_mag);

if met.valid
    fprintf('\n✓ GA ha trovato un flyby fisicamente valido.\n');
else
    fprintf('\n✗ GA non ha chiuso perfettamente il flyby. Serve refinement su TOF a path fissato.\n');
end
