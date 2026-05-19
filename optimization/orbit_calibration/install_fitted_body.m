function install_fitted_body(target_id, bodies_fit_file)
% Install a fitted body into the baseline bodies.mat with a backup.

if nargin < 2
    error('Usage: install_fitted_body(target_id, bodies_fit_file)');
end

this_dir = fileparts(mfilename('fullpath'));
opt_dir = fileparts(this_dir);
if ~exist('SETUP_mythara_paths_auto', 'file')
    addpath(opt_dir, '-begin');
end
cfg = SETUP_mythara_paths_auto();
backup_dir = cfg.backup_dir;
bodies_file = cfg.bodies_file;
project_root = fileparts(cfg.root);
data_bodies_file = fullfile(project_root, 'data', 'bodies.mat');
flyby_bodies_file = fullfile(cfg.flyby_dir, 'bodies.mat');

if ~isfile(bodies_file)
    error('bodies.mat not found: %s', bodies_file);
end
if ~isfile(bodies_fit_file)
    candidate = fullfile(cfg.calib_dir, 'bodies_fit', bodies_fit_file);
    if isfile(candidate)
        bodies_fit_file = candidate;
    else
        error('bodies_fit_file not found: %s', bodies_fit_file);
    end
end

load(bodies_file, 'bodies');
load(bodies_fit_file, 'bodies_fit');

key = sprintf('id_%d', target_id);
if ~isfield(bodies, key)
    error('Target ID %d not found in bodies.mat.', target_id);
end
if ~isfield(bodies_fit, key)
    error('Target ID %d not found in bodies_fit.', target_id);
end

if ~isfolder(backup_dir)
    mkdir(backup_dir);
end

backup_file = fullfile(backup_dir, sprintf('bodies_backup_%s.mat', datestr(now, 'yyyymmdd_HHMMSS')));
save(backup_file, 'bodies');

bodies.(key) = bodies_fit.(key);
save(bodies_file, 'bodies');

fprintf('Backup saved: %s\n', backup_file);
fprintf('Updated bodies.mat: %s\n', bodies_file);

extra_files = {data_bodies_file, flyby_bodies_file};
for k = 1:numel(extra_files)
    if isfile(extra_files{k})
        save(extra_files{k}, 'bodies');
        fprintf('Updated bodies.mat: %s\n', extra_files{k});
    else
        warning('Optional bodies.mat not found: %s', extra_files{k});
    end
end

end
