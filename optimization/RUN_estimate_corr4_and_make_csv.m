% Change only the INPUT block below.
clear; clc

cfg = SETUP_mythara_paths_auto();

global BODY_CORR_OVERRIDE

%% ================= INPUT =================
target_id = 8;
source_id = 7;
t_target  = 80;
shift_km  = 1.0;
dir       = 0;

E0  = NaN;   % base error
ExP = NaN;   % error +X
EyP = NaN;   % error +Y
EzP = NaN;   % error +Z
%% =========================================

if any(isnan([E0 ExP EyP EzP]))
    error('Inserisci E0, ExP, EyP, EzP dal validatore.')
end

s = shift_km;

corr = [
    (E0^2 + s^2 - ExP^2)/(2*s), ...
    (E0^2 + s^2 - EyP^2)/(2*s), ...
    (E0^2 + s^2 - EzP^2)/(2*s)
];

fprintf('\nEstimated correction:\n')
fprintf('id_%d, t = %.6f yr\n', target_id, t_target)
fprintf('C = [%.9f %.9f %.9f] km\n', corr)

solution.sequence    = [source_id target_id];
solution.t_flyby     = [0 t_target];
solution.Lambert_dir = dir;

BODY_CORR_OVERRIDE = struct();
BODY_CORR_OVERRIDE.(sprintf('id_%d',target_id)) = corr;

tag = time_tag(t_target);
outfile = fullfile(cfg.output_probe_dir, ...
    sprintf('corr_id%d_t%s.csv',target_id,tag));

solution_to_csv_patched(solution,outfile);

BODY_CORR_OVERRIDE = [];

fprintf('\nUpload corrected CSV:\n%s\n',outfile)

%% =====================================================================
function tag = time_tag(t)
tag = sprintf('%.6f',t);
tag = regexprep(tag,'0+$','');
tag = regexprep(tag,'\.$','');
tag = strrep(tag,'.','p');
end
