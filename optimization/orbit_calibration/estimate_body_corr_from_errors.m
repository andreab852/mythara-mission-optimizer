function [corr, outfile] = estimate_body_corr_from_errors(target_id, t_target, source_id, shift_km, E0, Ex_plus, Ex_minus, Ey_plus, Ey_minus, Ez_plus, Ez_minus)
% Estimate correction vector from scalar errors and write a corrected probe CSV.

if nargin < 11
    error(['Usage: estimate_body_corr_from_errors(target_id, t_target, source_id, shift_km, ', ...
        'E0, Ex_plus, Ex_minus, Ey_plus, Ey_minus, Ez_plus, Ez_minus)']);
end
if shift_km <= 0
    error('shift_km must be > 0.');
end

this_dir = fileparts(mfilename('fullpath'));
opt_dir = fileparts(this_dir);
if ~exist('SETUP_mythara_paths_auto', 'file')
    addpath(opt_dir, '-begin');
end
cfg = SETUP_mythara_paths_auto();

cx = (Ex_minus^2 - Ex_plus^2) / (4 * shift_km);
cy = (Ey_minus^2 - Ey_plus^2) / (4 * shift_km);
cz = (Ez_minus^2 - Ez_plus^2) / (4 * shift_km);

corr = [cx cy cz];

fprintf('E0 = %.6f km\n', E0);
fprintf('Estimated correction = [% .6f % .6f % .6f] km\n', corr(1), corr(2), corr(3));

out_dir = cfg.output_probe_dir;

solution.sequence = [source_id target_id];
solution.t_flyby = [0 t_target];
solution.Lambert_dir = 1;

tag = format_time_tag(t_target);
outfile = fullfile(out_dir, sprintf('corr_id%d_t%s.csv', target_id, tag));

apply_body_override(target_id, corr);

try
    evalc("solution_to_csv_patched(solution, outfile);");
catch ME
    clear_body_override();
    rethrow(ME);
end

clear_body_override();

fprintf('Corrected probe written: %s\n', outfile);

end

function apply_body_override(target_id, corr)

global BODY_CORR_OVERRIDE

BODY_CORR_OVERRIDE = struct();
BODY_CORR_OVERRIDE.(sprintf('id_%d', target_id)) = corr(:).';

end

function clear_body_override()

global BODY_CORR_OVERRIDE
BODY_CORR_OVERRIDE = [];

end

function tag = format_time_tag(t_year)

if abs(t_year - round(t_year)) < 1e-9
    tag = sprintf('%d', round(t_year));
else
    tag = sprintf('%.6f', t_year);
    tag = regexprep(tag, '0+$', '');
    tag = regexprep(tag, '\.$', '');
    tag = strrep(tag, '.', 'p');
end

end
