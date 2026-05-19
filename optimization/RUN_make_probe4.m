% Change only the INPUT block below.
clear; clc

cfg = SETUP_mythara_paths_auto();

%% ================= INPUT =================
target_id = 8;       % planet to calibrate
source_id = 7;       % planet already corrected
t_target  = 80;      % years
shift_km  = 1.0;     % km
%% =========================================

make_body_probe_4(target_id,t_target,source_id,shift_km);

fprintf('\nGenerated 4 probe files.\n')
fprintf('Target id_%d at t = %.6f yr\n', target_id, t_target)
fprintf('Upload: base, px, py, pz\n')
fprintf('Folder:\n%s\n', cfg.output_probe_dir)
