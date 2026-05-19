function bodies = generate_bodies(csv_file, out_file)
%GENERATE_BODIES Build Mythara bodies.mat from CSV/XLSX ephemeris.
%
% Usage:
%   bodies = generate_bodies();
%   bodies = generate_bodies('MytharaSystemData.csv');
%   bodies = generate_bodies('MytharaSystemData.csv','bodies.mat');
%
% This creates a struct compatible with solution_to_csv.m:
%   bodies.mu_star
%   bodies.id_4, bodies.id_5, ..., bodies.id_1000

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(csv_file)
    c1 = fullfile(here, 'MytharaSystemData.csv');
    c2 = fullfile(here, 'MytharaSystemData.xlsx');
    if isfile(c1)
        csv_file = c1;
    elseif isfile(c2)
        csv_file = c2;
    else
        error('Could not find MytharaSystemData.csv or MytharaSystemData.xlsx next to generate_bodies.m');
    end
end
if nargin < 2 || isempty(out_file)
    out_file = fullfile(here, 'bodies.mat');
end

T = local_read_ephemeris_table(csv_file);

required = {'ID','Name','GM','Radius','SMA','Ecc','Inc','RAAN','ArgP','M0','Weight'};
for k = 1:numel(required)
    if ~ismember(required{k}, T.Properties.VariableNames)
        error('Missing required column after import: %s', required{k});
    end
end

bodies = struct();
bodies.AU        = 149597870.691;     % km
bodies.mu_star   = 139348062043.343;  % km^3/s^2
bodies.day_sec   = 86400;
bodies.year_sec  = 365.25 * 86400;
bodies.valid_ids = [4 5 6 7 8 9 10 1000];

for i = 1:height(T)
    id = T.ID(i);
    field = sprintf('id_%d', id);

    body = struct();
    body.id     = id;
    body.name   = char(T.Name(i));

    body.a      = T.SMA(i);
    body.e      = T.Ecc(i);
    body.i      = T.Inc(i);

    body.omega  = T.RAAN(i);
    body.raan   = T.RAAN(i);

    body.w      = T.ArgP(i);
    body.argp   = T.ArgP(i);

    body.M0     = T.M0(i);

    body.gm     = T.GM(i);
    body.mu     = T.GM(i);

    body.radius = T.Radius(i);
    body.weight = T.Weight(i);

    body.i_rad     = deg2rad(body.i);
    body.omega_rad = deg2rad(body.omega);
    body.raan_rad  = deg2rad(body.raan);
    body.w_rad     = deg2rad(body.w);
    body.argp_rad  = deg2rad(body.argp);
    body.M0_rad    = deg2rad(body.M0);

    body.n_rad_s     = sqrt(bodies.mu_star / body.a^3);
    body.period_sec  = 2*pi / body.n_rad_s;
    body.period_year = body.period_sec / bodies.year_sec;

    bodies.(field) = body;
end

save(out_file, 'bodies');
fprintf('Saved %d bodies to %s\n', height(T), out_file);
end

function T = local_read_ephemeris_table(file)
[~,~,ext] = fileparts(file);
if strcmpi(ext,'.csv')
    raw = readtable(file, 'Delimiter',';', 'VariableNamingRule','preserve', 'TextType','string');
elseif any(strcmpi(ext,{'.xlsx','.xls'}))
    raw = readtable(file, 'VariableNamingRule','preserve', 'TextType','string');
else
    error('Unsupported ephemeris file type: %s', ext);
end

% Rename columns explicitly from the actual Mythara header names
orig = string(raw.Properties.VariableNames);
ren  = orig;
for i = 1:numel(orig)
    s = strtrim(orig(i));
    s = erase(s, '#');
    switch s
        case 'Planet ID'
            ren(i) = "ID";
        case 'Name'
            ren(i) = "Name";
        case 'GM (km3/s2)'
            ren(i) = "GM";
        case 'Radius (km)'
            ren(i) = "Radius";
        case 'Semi-Major Axis (km)'
            ren(i) = "SMA";
        case 'Eccentricity ()'
            ren(i) = "Ecc";
        case 'Inclination (deg)'
            ren(i) = "Inc";
        case 'Longitude of the Ascending Node (deg)'
            ren(i) = "RAAN";
        case 'Argument of Periapsis (deg)'
            ren(i) = "ArgP";
        case 'Mean Anomaly at t=0 (deg)'
            ren(i) = "M0";
        case 'Weight ()'
            ren(i) = "Weight";
    end
end
raw.Properties.VariableNames = cellstr(ren);

want = {'ID','Name','GM','Radius','SMA','Ecc','Inc','RAAN','ArgP','M0','Weight'};
T = raw(:, want);

% Convert numeric columns
numvars = setdiff(want, {'Name'});
for i = 1:numel(numvars)
    v = numvars{i};
    T.(v) = str2double(string(T.(v)));
end
T.Name = string(T.Name);

% Remove empty/spacer rows
T = T(~isnan(T.ID), :);
end
