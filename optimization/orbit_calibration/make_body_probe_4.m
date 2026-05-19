function make_body_probe_4(target_id,t_target,source_id,shift_km)

global BODY_CORR_OVERRIDE

if nargin < 4
    shift_km = 1.0;
end
this_dir = fileparts(mfilename('fullpath'));
opt_dir = fileparts(this_dir);
if ~exist('SETUP_mythara_paths_auto', 'file')
    addpath(opt_dir, '-begin');
end
cfg = SETUP_mythara_paths_auto();
outdir = cfg.output_probe_dir;

tags = {'base','px','py','pz'};
corrs = [
    0 0 0
    shift_km 0 0
    0 shift_km 0
    0 0 shift_km
];

for k = 1:numel(tags)

    solution.sequence    = [source_id target_id];
    solution.t_flyby     = [0 t_target];
    solution.Lambert_dir = 0;

    BODY_CORR_OVERRIDE = struct();
    BODY_CORR_OVERRIDE.(sprintf('id_%d',target_id)) = corrs(k,:);

    tag = format_time_tag(t_target);
    outfile = fullfile(outdir, ...
        sprintf('probe_id%d_t%s_%s.csv',target_id,tag,tags{k}));

    solution_to_csv_patched(solution,outfile);

    fprintf('%s | corr=[%.6f %.6f %.6f]\n', ...
        outfile,corrs(k,1),corrs(k,2),corrs(k,3))
end

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
