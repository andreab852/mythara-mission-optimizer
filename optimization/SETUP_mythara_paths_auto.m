function cfg = SETUP_mythara_paths_auto()
%SETUP_MYTHARA_PATHS_AUTO Canonical setup for new optimization folder layout.

root = fileparts(mfilename('fullpath'));

cfg = struct();
cfg.root = root;

cfg.flyby_dir = fullfile(root,'flyby_functions');
cfg.calib_dir = fullfile(root,'orbit_calibration');
cfg.submission_dir = fullfile(root,'submission');
cfg.utils_dir = fullfile(cfg.submission_dir,'utils');
cfg.backup_dir = fullfile(cfg.submission_dir,'Backup');
cfg.search_dir = cfg.submission_dir;

cfg.KM_solver = fullfile(cfg.flyby_dir,'KM_solver.m');
cfg.lambert_universal = fullfile(cfg.flyby_dir,'lambert_universal.m');
cfg.stumpff = fullfile(cfg.flyby_dir,'stumpff.m');

cfg.solution_to_csv_patched = fullfile(cfg.utils_dir,'solution_to_csv_patched.m');
cfg.bodies_file = fullfile(cfg.submission_dir,'bodies.mat');

cfg.output_probe_dir = fullfile(cfg.calib_dir,'output_probe');

required = {
    cfg.flyby_dir
    cfg.calib_dir
    cfg.submission_dir
    cfg.utils_dir
    cfg.KM_solver
    cfg.lambert_universal
    cfg.stumpff
    cfg.solution_to_csv_patched
    cfg.bodies_file
};

for k = 1:numel(required)
    if ~(isfolder(required{k}) || isfile(required{k}))
        error('Missing required path/file:\n%s', required{k})
    end
end

addpath(cfg.flyby_dir,'-begin')
addpath(cfg.calib_dir,'-begin')
addpath(cfg.submission_dir,'-begin')
addpath(cfg.utils_dir,'-begin')
addpath(root,'-begin')

if ~exist(cfg.output_probe_dir,'dir')
    mkdir(cfg.output_probe_dir)
end
if ~exist(cfg.backup_dir,'dir')
    mkdir(cfg.backup_dir)
end

rehash toolboxcache

fprintf('\n===== MYTHARA SETUP OK =====\n')
fprintf('root:       %s\n',cfg.root)
fprintf('flyby:      %s\n',cfg.flyby_dir)
fprintf('calib:      %s\n',cfg.calib_dir)
fprintf('submit:     %s\n',cfg.submission_dir)
fprintf('utils:      %s\n',cfg.utils_dir)
fprintf('bodies:     %s\n',cfg.bodies_file)

fprintf('\nUsing:\n')
which KM_solver
which lambert_universal
which stumpff
which solution_to_csv_patched
which make_body_probe_4
which fit_body_orbit_from_points
which install_fitted_body

end