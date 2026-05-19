% Change only the INPUT block below.
clear; clc

cfg = SETUP_mythara_paths_auto();

%% ================= INPUT =================
target_id = 8;
model = 'full';

T = [80; 160; 240];

C = [
   -1.729490  -0.178130   0.217716
    2.846213  -1.126220  -0.674329
   -1.716605   4.100281   1.168870
];
%% =========================================

[bodies_fit,bfit,rms] = fit_body_orbit_from_points(target_id,T,C,model);

fprintf('\nFit complete for id_%d\n',target_id)
fprintf('RMS = %.9f km\n',rms)

fprintf('\nIf RMS < 0.1 km, install with:\n')
fprintf('install_fitted_body(%d,''bodies_fit_id_%d.mat'')\n',target_id,target_id)
fprintf('\nFit file should be inside orbit_calibration folder.\n')
fprintf('Setup bodies used:\n%s\n',cfg.bodies_file)
