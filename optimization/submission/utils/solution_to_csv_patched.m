function solution_to_csv_patched(solution, output_file)
% Patched copy of solution_to_csv.
% Corrects selected planet states before Lambert and CSV writing.

solution = normalize_input(solution);

seq = solution.sequence(:)';
n_planets = length(seq);
n_legs = n_planets - 1;
t_flyby = solution.t_flyby(:)';
Lambert_dir = round(solution.Lambert_dir(:)');

this_dir = fileparts(mfilename('fullpath'));
submission_dir = fileparts(this_dir);
bodies_file = fullfile(submission_dir, 'bodies.mat');

if ~isfile(bodies_file)
    bodies_file = fullfile(this_dir, 'bodies.mat');
end
if ~isfile(bodies_file)
    error('bodies.mat not found in submission folder: %s', submission_dir);
end

tmp = load(bodies_file, 'bodies');
bodies = tmp.bodies;

mu = bodies.mu_star;
if isfield(bodies,'year_sec')
    y2sec = bodies.year_sec;
else
    y2sec = 365.25*24*3600;
end

for k = 1:n_planets
    field_name = sprintf('id_%d', seq(k));
    if ~isfield(bodies, field_name)
        error('Planet ID %d not found in bodies.mat.', seq(k));
    end
end

r_planets = zeros(n_planets,3);
v_planets = zeros(n_planets,3);

for k = 1:n_planets
    state_corr = body_state_patched(seq(k), t_flyby(k), bodies, y2sec);
    r_planets(k,:) = state_corr(1:3);
    v_planets(k,:) = state_corr(4:6);
end

v_dep = zeros(n_legs,3);
v_arr = zeros(n_legs,3);

for leg = 1:n_legs
    tof_sec = (t_flyby(leg+1) - t_flyby(leg))*y2sec;

    if tof_sec <= 0
        error('Non-positive TOF for leg %d.', leg);
    end

    r1 = r_planets(leg,:)';
    r2 = r_planets(leg+1,:)';

    tA = t_flyby(leg)*y2sec;
    tB = t_flyby(leg+1)*y2sec;

    try
        [vd, va, ~] = lambert_universal(mu, tA, r1, tB, r2, Lambert_dir(leg));

        % Final differential correction on departure velocity.
        % This forces the conic propagated from row A to hit the exact row-B
        % position written in the CSV. Needed because the validator tolerance
        % is 0.1 km and small Lambert residuals of 0.14-0.40 km are enough to fail.
        [vd, va, pos_err, dv_fix] = correct_lambert_arc_dc(r1, r2, real(vd(:)), tof_sec, mu);

        v_dep(leg,:) = real(vd(:))';
        v_arr(leg,:) = real(va(:))';

        if pos_err > 0.05
            fprintf('WARNING: Leg %d (%d->%d) position error = %.6f km (TOF=%.2f yr, |dv_fix|=%.3e km/s)\n', ...
                leg, seq(leg), seq(leg+1), pos_err, tof_sec/y2sec, dv_fix);
        end

    catch e
        error('Lambert solver failed for leg %d (%d -> %d, TOF=%.2f yr): %s', ...
            leg, seq(leg), seq(leg+1), t_flyby(leg+1)-t_flyby(leg), e.message);
    end
end

v_inf_in  = zeros(n_planets,3);
v_inf_out = zeros(n_planets,3);

for k = 1:n_planets
    v_p = v_planets(k,:);

    if k == 1
        if n_legs > 0
            v_inf_in(k,:) = v_dep(1,:) - v_p;
        end
    else
        v_inf_in(k,:) = v_arr(k-1,:) - v_p;
    end

    if k <= n_legs
        v_inf_out(k,:) = v_dep(k,:) - v_p;
    end
end

fid = fopen(output_file, 'w');
if fid == -1
    error('Cannot open file: %s', output_file);
end

fprintf(fid, '# Patched GTOC13 Solution File\n');
fprintf(fid, '# Generated: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, '# Sequence: [%s]\n', sprintf('%d ', seq));
fprintf(fid, '# Mission time: %.2f years\n', t_flyby(end) - t_flyby(1));
fprintf(fid, '#\n');

for k = 1:n_planets
    epoch_sec = t_flyby(k)*y2sec;
    r = r_planets(k,:);

    if k <= n_legs
        v_sc_in  = v_planets(k,:) + v_inf_in(k,:);
        v_sc_out = v_dep(k,:);

        write_row(fid, seq(k), 1, epoch_sec, r, v_sc_in,  v_inf_in(k,:));
        write_row(fid, seq(k), 1, epoch_sec, r, v_sc_out, v_inf_out(k,:));

        epoch_end = t_flyby(k+1)*y2sec;
        r_next = r_planets(k+1,:);

        write_row(fid, 0, 0, epoch_sec,  r,      v_dep(k,:), [0 0 0]);
        write_row(fid, 0, 0, epoch_end,  r_next, v_arr(k,:), [0 0 0]);

    else
        v_sc_in = v_planets(k,:) + v_inf_in(k,:);
        write_row(fid, seq(k), 1, epoch_sec, r, v_sc_in, v_inf_in(k,:));
    end
end

fclose(fid);

fprintf('Patched CSV written: %s\n', output_file);
fprintf('Sequence: [%s]\n', sprintf('%d ', seq));
fprintf('Mission time: %.2f years\n', t_flyby(end) - t_flyby(1));

try
    flyby_data = build_flyby_data(solution, r_planets, v_inf_in, v_inf_out, t_flyby);
    score_data = compute_gtoc13_score(flyby_data);
    fprintf('Score: %.4f\n', score_data.score);
    fprintf('Score/year: %.4f\n', score_data.score_per_year);
catch e
    fprintf('Warning: could not compute score: %s\n', e.message);
end

end

function state_corr = body_state_patched(id, t_year, bodies, y2sec)

state0 = KM_solver(t_year*y2sec, bodies.(sprintf('id_%d',id)));
r0 = state0(1:3).';
corr0 = body_corr(id, t_year, bodies);
r_corr = r0 + corr0;

dt = 10; % seconds

if t_year*y2sec > dt
    tp = t_year + dt/y2sec;
    tm = t_year - dt/y2sec;

    sp = KM_solver(tp*y2sec, bodies.(sprintf('id_%d',id)));
    sm = KM_solver(tm*y2sec, bodies.(sprintf('id_%d',id)));

    rp = sp(1:3).' + body_corr(id, tp, bodies);
    rm = sm(1:3).' + body_corr(id, tm, bodies);

    v_corr = (rp - rm)/(2*dt);
else
    tp = t_year + dt/y2sec;

    sp = KM_solver(tp*y2sec, bodies.(sprintf('id_%d',id)));

    rp = sp(1:3).' + body_corr(id, tp, bodies);

    v_corr = (rp - r_corr)/dt;
end

state_corr = [r_corr v_corr];

end

function corr = body_corr(id, t_year, bodies)

global BODY_CORR_OVERRIDE

key = sprintf('id_%d',id);

if ~isempty(BODY_CORR_OVERRIDE) && ...
        isstruct(BODY_CORR_OVERRIDE) && ...
        isfield(BODY_CORR_OVERRIDE,key)

    corr = BODY_CORR_OVERRIDE.(key);
    return
end

corr = [0 0 0];

switch id
    case 5
        corr = [0 0 0];

    case 4
        corr = [0 0 0];
end

end


function [v1, v2, pos_err, dv_total] = correct_lambert_arc_dc(r1, r2, v1, dt, mu)
% Robust Newton/shooting correction on v1 so the exported conic arc lands
% on the exact exported arrival position. This is only an export-level fix.
% Typical corrections are << 0.1 mm/s for 0.1-1 km residuals over years.

r1 = r1(:);
r2 = r2(:);
v1 = v1(:);

dv_total = 0;

[r_end, v_end] = kepler_propagate_fg_state(r1, v1, dt, mu);
pos_err = norm(r_end - r2);
if ~isfinite(pos_err)
    v2 = v_end;
    return
end

% Do not silently repair kilometer-to-thousands-km arcs: those are bad
% solutions / bad ephemeris corrections, not Lambert roundoff.
if pos_err > 5.0
    v2 = v_end;
    return
end

% Adaptive shooting. The previous version used a too-small Jacobian step and
% an overly tight dv cap, so 0.18-0.37 km residuals could survive unchanged.
for it = 1:12
    [r_end, ~] = kepler_propagate_fg_state(r1, v1, dt, mu);
    err = r_end - r2;
    pos_err = norm(err);

    if pos_err < 1e-4      % 0.1 m, validator tolerance is 100 m
        break
    end

    % Finite-difference scale: not too small, otherwise universal-variable
    % propagation noise makes J singular.
    dv_step = max(1e-7, 1e-8*max(1,norm(v1))); % km/s
    J = zeros(3,3);
    for j = 1:3
        e = zeros(3,1);
        e(j) = dv_step;
        rp = kepler_propagate_fg_state(r1, v1 + e, dt, mu);
        rm = kepler_propagate_fg_state(r1, v1 - e, dt, mu);
        J(:,j) = (rp - rm)/(2*dv_step);
    end

    if any(~isfinite(J),'all') || rcond(J) < 1e-16
        break
    end

    dv = -J \ err;
    if any(~isfinite(dv))
        break
    end

    % Allow enough correction for sub-km residuals over multi-year arcs.
    % 2e-6 km/s = 2 mm/s, still tiny dynamically; line search usually uses far less.
    max_dv_iter = 2e-6;
    ndv = norm(dv);
    if ndv > max_dv_iter
        dv = dv * (max_dv_iter/ndv);
    end

    % Backtracking: only accept if residual decreases.
    base_err = pos_err;
    accepted = false;
    alpha = 1.0;
    for ls = 1:10
        v_try = v1 + alpha*dv;
        r_try = kepler_propagate_fg_state(r1, v_try, dt, mu);
        try_err = norm(r_try - r2);
        if isfinite(try_err) && try_err < base_err
            v1 = v_try;
            dv_total = dv_total + norm(alpha*dv);
            accepted = true;
            break
        end
        alpha = 0.5*alpha;
    end

    if ~accepted
        break
    end

    if norm(alpha*dv) < 1e-13
        break
    end
end

[r_end, v2] = kepler_propagate_fg_state(r1, v1, dt, mu);
pos_err = norm(r_end - r2);

end

function s = normalize_input(s)

has_t_flyby = isfield(s, 't_flyby');
has_t1_TOF  = isfield(s, 't1') && isfield(s, 'TOF');

if ~isfield(s, 'sequence')
    error('Missing field: solution.sequence');
end
if ~isfield(s, 'Lambert_dir')
    error('Missing field: solution.Lambert_dir');
end

seq = s.sequence(:)';
n = length(seq);

if n < 2
    error('sequence must contain at least 2 planets.');
end

valid_ids = [4, 5, 6, 7, 8, 9, 10, 1000];
for k = 1:n
    if ~ismember(seq(k), valid_ids)
        error('Invalid planet ID %d.', seq(k));
    end
end

n_legs = n - 1;

if length(s.Lambert_dir) ~= n_legs
    error('Lambert_dir must have %d elements.', n_legs);
end

if any(s.Lambert_dir ~= 0 & s.Lambert_dir ~= 1)
    error('Lambert_dir values must be 0 or 1.');
end

if has_t_flyby
    if length(s.t_flyby) ~= n
        error('t_flyby must have %d elements.', n);
    end

    t_flyby = s.t_flyby(:)';

    if any(diff(t_flyby) <= 0)
        error('t_flyby must be strictly increasing.');
    end

elseif has_t1_TOF
    if ~isscalar(s.t1) || s.t1 < 0
        error('t1 must be a non-negative scalar.');
    end

    if length(s.TOF) ~= n_legs
        error('TOF must have %d elements.', n_legs);
    end

    if any(s.TOF <= 0)
        error('All TOF values must be positive.');
    end

    t_flyby = zeros(1,n);
    t_flyby(1) = s.t1;

    for k = 2:n
        t_flyby(k) = t_flyby(k-1) + s.TOF(k-1);
    end

    s.t_flyby = t_flyby;
else
    error('Must provide either t_flyby or t1 + TOF.');
end

s.sequence = seq;
s.t_flyby = t_flyby;

end

function fd = build_flyby_data(solution, r_planets, v_inf_in, v_inf_out, t_flyby)

fd.sequence = solution.sequence(:)';
fd.times = t_flyby;
fd.positions = r_planets;
fd.velocities_in_planet_frame  = v_inf_in;
fd.velocities_out_planet_frame = v_inf_out;
fd.v_inf_in  = vecnorm(v_inf_in, 2, 2)';
fd.v_inf_out = vecnorm(v_inf_out, 2, 2)';

end

function write_row(fid, body_id, flag, epoch, r, v, ctrl)

r = r(:)';
v = v(:)';
ctrl = ctrl(:)';

fprintf(fid, '%d, %d, %.15e, %.15e, %.15e, %.15e, %.15e, %.15e, %.15e, %.15e, %.15e, %.15e\n', ...
    body_id, flag, epoch, ...
    r(1), r(2), r(3), ...
    v(1), v(2), v(3), ...
    ctrl(1), ctrl(2), ctrl(3));

end

function r2 = kepler_propagate_fg(r1_vec, v1_vec, dt, mu)

r1_vec = r1_vec(:);
v1_vec = v1_vec(:);

r1 = norm(r1_vec);
v1 = norm(v1_vec);
vr1 = dot(r1_vec, v1_vec)/r1;
alpha = 2/r1 - v1^2/mu;

if alpha > 1e-10
    chi0 = sqrt(mu)*dt*alpha;
elseif alpha < -1e-10
    a = 1/alpha;
    chi0 = sign(dt)*sqrt(-a)*log(-2*mu*alpha*dt / ...
        (dot(r1_vec,v1_vec) + sign(dt)*sqrt(-mu*a)*(1-r1*alpha)));
else
    h = cross(r1_vec, v1_vec);
    p = norm(h)^2/mu;
    s = 0.5*atan(1/(3*sqrt(mu/p^3)*dt));
    w = atan(tan(s)^(1/3));
    chi0 = sqrt(p)*2/tan(2*w);
end

chi = chi0;

for k = 1:300
    psi = chi^2*alpha;
    [C, S] = stumpff(psi);

    r = chi^2*C + vr1/sqrt(mu)*chi*(1-psi*S) + r1*(1-psi*C);

    F = r1*vr1/sqrt(mu)*chi^2*C + ...
        (1-r1*alpha)*chi^3*S + r1*chi - sqrt(mu)*dt;

    chi = chi - F/r;

    if abs(F) < 1e-12*sqrt(mu)*abs(dt)
        break;
    end
end

psi = chi^2*alpha;
[C, S] = stumpff(psi);

f = 1 - chi^2/r1*C;
g = dt - chi^3/sqrt(mu)*S;

r2 = f*r1_vec + g*v1_vec;

end


function [r2, v2] = kepler_propagate_fg_state(r1_vec, v1_vec, dt, mu)

r1_vec = r1_vec(:);
v1_vec = v1_vec(:);

r1 = norm(r1_vec);
v1 = norm(v1_vec);
vr1 = dot(r1_vec, v1_vec)/r1;
alpha = 2/r1 - v1^2/mu;

if alpha > 1e-10
    chi0 = sqrt(mu)*dt*alpha;
elseif alpha < -1e-10
    a = 1/alpha;
    arg = -2*mu*alpha*dt / ...
        (dot(r1_vec,v1_vec) + sign(dt)*sqrt(-mu*a)*(1-r1*alpha));
    if arg <= 0 || ~isfinite(arg)
        chi0 = sign(dt)*sqrt(abs(a));
    else
        chi0 = sign(dt)*sqrt(-a)*log(arg);
    end
else
    h = cross(r1_vec, v1_vec);
    p = norm(h)^2/mu;
    s = 0.5*atan(1/(3*sqrt(mu/p^3)*dt));
    w = atan(tan(s)^(1/3));
    chi0 = sqrt(p)*2/tan(2*w);
end

chi = chi0;

for k = 1:300
    psi = chi^2*alpha;
    [C, S] = stumpff(psi);

    r = chi^2*C + vr1/sqrt(mu)*chi*(1-psi*S) + r1*(1-psi*C);

    F = r1*vr1/sqrt(mu)*chi^2*C + ...
        (1-r1*alpha)*chi^3*S + r1*chi - sqrt(mu)*dt;

    if ~isfinite(F) || ~isfinite(r) || abs(r) < eps
        break
    end

    chi = chi - F/r;

    if abs(F) < 1e-13*sqrt(mu)*abs(dt)
        break;
    end
end

psi = chi^2*alpha;
[C, S] = stumpff(psi);

f = 1 - chi^2/r1*C;
g = dt - chi^3/sqrt(mu)*S;

r2 = f*r1_vec + g*v1_vec;
r2n = norm(r2);

fdot = sqrt(mu)/(r2n*r1) * chi * (psi*S - 1);
gdot = 1 - chi^2/r2n*C;

v2 = fdot*r1_vec + gdot*v1_vec;

end
