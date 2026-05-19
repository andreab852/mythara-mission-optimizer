function files = make_body_probe_batch(target_id, t_target, source_id, shift_km)
% Generate 7 probe CSVs with +/- axis shifts on the target body.

if nargin < 3
    error('Usage: make_body_probe_batch(target_id, t_target, source_id, shift_km)');
end
if nargin < 4 || isempty(shift_km)
    shift_km = 1.0;
end
if t_target <= 0
    error('t_target must be > 0 years.');
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
out_dir = cfg.output_probe_dir;

solution.sequence = [source_id target_id];
solution.t_flyby = [0 t_target];
solution.Lambert_dir = 1;

tag = format_time_tag(t_target);
labels = {'base','px','mx','py','my','pz','mz'};
shifts = [
    0 0 0
    1 0 0
   -1 0 0
    0 1 0
    0 -1 0
    0 0 1
    0 0 -1
];

files = cell(1, numel(labels));

for k = 1:numel(labels)
    corr = shift_km * shifts(k,:);
    outfile = fullfile(out_dir, sprintf('probe_id%d_t%s_%s.csv', target_id, tag, labels{k}));

    apply_body_override(target_id, corr);
    evalc("solution_to_csv_patched(solution, outfile);");

    files{k} = outfile;
end

clear_body_override();

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
