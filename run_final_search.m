%RUN_FINAL_SEARCH Launch the final Mythara multi-flyby search workflow.
% This wrapper is robust to the current MATLAB folder: it resolves paths
% from its own location and then runs the final parallel script.

clear; clc

project_root = fileparts(mfilename('fullpath'));
script_file = fullfile(project_root, 'optimization', 'submission', ...
    'Main_Flyby_multistrategy.m');

run(script_file);
