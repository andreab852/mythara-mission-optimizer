function out = compute_gtoc13_score(flyby_data)
%COMPUTE_GTOC13_SCORE Mythara project scoring helper.
%
% Expects a struct with fields like those built by solution_to_csv:
%   sequence     1xN body IDs
%   times        1xN flyby epochs [years] (optional for score_per_year)
%   positions    Nx3 heliocentric positions [km]
%   v_inf_in     1xN incoming hyperbolic excess speeds [km/s]
%
% Optional:
%   scientific_flags 1xN logical/double. If absent, all flybys are scientific.
%
% Implemented from the project brief:
%   J = b * sum_k sum_i w_k * S_{k,i} * F(Vinf_{k,i})
%   b = 1.2 if all bodies {4,5,6,7,8,9,10,1000} are visited scientifically,
%       else 1.0
%   S_1 = 1
%   S_i = 0.1 + 0.9 / (1 + 10 * sum_{j<i} exp(-(acosd(r_i·r_j)/50)^2))
%   F(v) = 0.2 + exp(-v/13)/(1 + exp(-5*(v-1.5)))
%
% Returns fields:
%   score, score_per_year, raw_sum, bonus, n_scientific, body_breakdown

weights = containers.Map('KeyType','double','ValueType','double');
weights(4) = 3;
weights(1000) = 5;
weights(5) = 7;
weights(6) = 10;
weights(7) = 15;
weights(8) = 20;
weights(9) = 35;
weights(10) = 50;
all_ids = [4 5 6 7 8 9 10 1000];

seq = flyby_data.sequence(:)';
N = numel(seq);
if isfield(flyby_data, 'scientific_flags')
    sci = logical(flyby_data.scientific_flags(:)');
else
    sci = true(1,N);
end
if numel(sci) ~= N
    error('scientific_flags must have the same length as sequence.');
end

if isfield(flyby_data, 'positions')
    R = flyby_data.positions;
else
    error('compute_gtoc13_score requires flyby_data.positions (Nx3).');
end
if size(R,1) ~= N || size(R,2) ~= 3
    error('positions must be Nx3.');
end

if isfield(flyby_data, 'v_inf_in')
    vinf = flyby_data.v_inf_in(:)';
elseif isfield(flyby_data, 'velocities_in_planet_frame')
    vinf = vecnorm(flyby_data.velocities_in_planet_frame, 2, 2)';
else
    error('compute_gtoc13_score requires v_inf_in or velocities_in_planet_frame.');
end
if numel(vinf) ~= N
    error('v_inf_in must have the same length as sequence.');
end

% Normalize heliocentric directions.
rhat = R ./ max(vecnorm(R, 2, 2), eps);

visited_science = unique(seq(sci), 'stable');
bonus = 1.0;
if all(ismember(all_ids, visited_science))
    bonus = 1.2;
end

raw_sum = 0;
body_breakdown = struct();
for id = all_ids
    idx = find(seq == id & sci);
    if isempty(idx)
        continue;
    end
    if numel(idx) > 13
        idx = idx(1:13);
    end

    contrib = zeros(1, numel(idx));
    prev_dirs = [];
    for ii = 1:numel(idx)
        k = idx(ii);
        if isempty(prev_dirs)
            S = 1.0;
        else
            dots = prev_dirs * rhat(k,:)';
            dots = max(-1, min(1, dots));
            ang = real(acosd(dots));
            S = 0.1 + 0.9/(1 + 10*sum(exp(-(ang.^2)/50)));
        end
        F = 0.2 + exp(-vinf(k)/13) / (1 + exp(-5*(vinf(k) - 1.5)));
        contrib(ii) = weights(id) * S * F;
        prev_dirs = [prev_dirs; rhat(k,:)]; %#ok<AGROW>
    end

    raw_sum = raw_sum + sum(contrib);
    body_breakdown.(sprintf('id_%d', id)) = struct( ...
        'count_scientific', numel(idx), ...
        'contribution_unscaled', sum(contrib), ...
        'contributions_each', contrib);
end

score = bonus * raw_sum;
if isfield(flyby_data, 'times') && numel(flyby_data.times) >= 2
    mission_years = flyby_data.times(end) - flyby_data.times(1);
else
    mission_years = NaN;
end

if isfinite(mission_years) && mission_years > 0
    score_per_year = score / mission_years;
else
    score_per_year = NaN;
end

out = struct();
out.score = score;
out.score_per_year = score_per_year;
out.raw_sum = raw_sum;
out.bonus = bonus;
out.n_scientific = sum(sci);
out.body_breakdown = body_breakdown;
end
