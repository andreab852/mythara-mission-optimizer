function [bodies_fit, bfit, rms] = fit_body_orbit_from_points(target_id, T, C, model)
% Fit orbital elements from corrected inertial points.

if nargin < 3
    error('Usage: fit_body_orbit_from_points(target_id, T, C, model)');
end
if nargin < 4 || isempty(model)
    model = 'full';
end

this_dir = fileparts(mfilename('fullpath'));
opt_dir = fileparts(this_dir);
if ~exist('SETUP_mythara_paths_auto', 'file')
    addpath(opt_dir, '-begin');
end
cfg = SETUP_mythara_paths_auto();

T = T(:);
if size(C,1) ~= numel(T) || size(C,2) ~= 3
    error('C must be an Nx3 matrix aligned with T.');
end
bodies_file = cfg.bodies_file;
load(bodies_file, 'bodies');

key = sprintf('id_%d', target_id);
if ~isfield(bodies, key)
    error('Target ID %d not found in bodies.mat.', target_id);
end

b0 = bodies.(key);
mu = bodies.mu_star;
if isfield(bodies, 'year_sec')
    ys = bodies.year_sec;
else
    ys = 365.25 * 24 * 3600;
end

fprintf('\n===== FIT BODY ORBIT (id_%d) =====\n', target_id);

if strcmpi(model, 'auto')
    [x1, rms1, b1] = fit_model('aM0', b0, T, C, ys, mu);
    [x2, rms2, b2] = fit_model('full', b0, T, C, ys, mu);
    if rms2 < rms1
        model = 'full';
        rms = rms2;
        bfit = b2;
    else
        model = 'aM0';
        rms = rms1;
        bfit = b1;
    end
    fprintf('\nBest model: %s (RMS %.9f km)\n', model, rms);
else
    [x, rms, bfit] = fit_model(model, b0, T, C, ys, mu); %#ok<ASGLU>
end

fprintf('\nDelta values:\n');
fprintf('da       = %.12f km\n', bfit.a - b0.a);
fprintf('de       = %.12e\n', bfit.e - b0.e);
fprintf('dargp    = %.12e deg\n', rad2deg(bfit.argp_rad - b0.argp_rad));
fprintf('dM0      = %.12e deg\n', rad2deg(bfit.M0_rad - b0.M0_rad));

fprintf('\nCheck residuals at points:\n');
for k = 1:numel(T)
    s0 = KM_solver(T(k) * ys, b0);
    sf = KM_solver(T(k) * ys, bfit);

    r_true = s0(1:3).' + C(k,:);
    err = norm(sf(1:3).' - r_true);

    fprintf('t = %8.3f yr | err = %.9f km\n', T(k), err);
end

bodies_fit = bodies;
bodies_fit.(key) = bfit;

out_dir = fullfile(cfg.calib_dir, 'bodies_fit');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
out_file = fullfile(out_dir, sprintf('bodies_fit_id_%d.mat', target_id));
save(out_file, 'bodies_fit');

fprintf('\nSaved: %s\n', out_file);

end

function [x, rms, bfit] = fit_model(mode, b0, T, C, ys, mu)

switch lower(mode)
    case 'am0'
        x = [0; 0];
        h = [0.1; 1.0];

    case 'full'
        x = [0; 0; 0; 0];
        h = [0.1; 1.0; 1.0; 1.0];

    otherwise
        error('Unknown model: %s', mode);
end

for iter = 1:30

    res0 = residual_vec(x, mode, b0, T, C, ys, mu);
    rms0 = sqrt(mean(res0.^2));

    J = finite_jacobian(x, mode, b0, T, C, ys, mu, h);

    JTJ = J.' * J;
    g   = J.' * res0;

    lam = 1e-10 * max(1, norm(JTJ, 'fro'));
    dx = -(JTJ + lam * eye(numel(x))) \ g;

    alpha = 1.0;
    accepted = false;

    for ls = 1:20
        xt = x + alpha * dx;

        if is_bad_params(xt, mode, b0)
            alpha = alpha / 2;
            continue
        end

        rest = residual_vec(xt, mode, b0, T, C, ys, mu);
        rmst = sqrt(mean(rest.^2));

        if rmst < rms0
            x = xt;
            accepted = true;
            break
        end

        alpha = alpha / 2;
    end

    fprintf('iter %02d | rms %.9f km | step %.3e | alpha %.3f\n', ...
        iter, rms0, norm(dx), alpha);

    if ~accepted || norm(alpha * dx) < 1e-10
        break
    end
end

res = residual_vec(x, mode, b0, T, C, ys, mu);
rms = sqrt(mean(res.^2));
bfit = apply_x(b0, x, mode, mu, ys);

fprintf('\nFinal %s RMS = %.9f km\n', mode, rms);

end

function J = finite_jacobian(x, mode, b0, T, C, ys, mu, h)

r0 = residual_vec(x, mode, b0, T, C, ys, mu); %#ok<NASGU>
J = zeros(3 * numel(T), numel(x));

for j = 1:numel(x)
    xp = x;
    xm = x;

    xp(j) = xp(j) + h(j);
    xm(j) = xm(j) - h(j);

    rp = residual_vec(xp, mode, b0, T, C, ys, mu);
    rm = residual_vec(xm, mode, b0, T, C, ys, mu);

    J(:, j) = (rp - rm) / (2 * h(j));
end

end

function res = residual_vec(x, mode, b0, T, C, ys, mu)

b = apply_x(b0, x, mode, mu, ys);
res = zeros(3 * numel(T), 1);

for k = 1:numel(T)
    s_nom = KM_solver(T(k) * ys, b0);
    s_fit = KM_solver(T(k) * ys, b);

    r_true = s_nom(1:3).' + C(k,:);
    r_fit  = s_fit(1:3).';

    idx = 3 * k - 2:3 * k;
    res(idx) = (r_fit - r_true).';
end

end

function b = apply_x(b0, x, mode, mu, ys)

b = b0;

switch lower(mode)
    case 'am0'
        da       = x(1);
        dM0_rad  = x(2) * 1e-9;

        b.a      = b0.a + da;
        b.M0_rad = b0.M0_rad + dM0_rad;

    case 'full'
        da        = x(1);
        de        = x(2) * 1e-9;
        dargp_rad = x(3) * 1e-9;
        dM0_rad   = x(4) * 1e-9;

        b.a        = b0.a + da;
        b.e        = b0.e + de;
        b.argp_rad = b0.argp_rad + dargp_rad;
        b.w_rad    = b.argp_rad;
        b.M0_rad   = b0.M0_rad + dM0_rad;
end

b.argp = rad2deg(b.argp_rad);
b.w    = b.argp;
b.M0   = rad2deg(b.M0_rad);

b.n_rad_s     = sqrt(mu / b.a^3);
b.period_sec  = 2 * pi / b.n_rad_s;
b.period_year = b.period_sec / ys;

end

function bad = is_bad_params(x, mode, b0)

bad = false;

switch lower(mode)
    case 'am0'
        a = b0.a + x(1);
        e = b0.e;

    case 'full'
        a = b0.a + x(1);
        e = b0.e + x(2) * 1e-9;

    otherwise
        error('Unknown model: %s', mode);
end

if ~isfinite(a) || ~isfinite(e) || a <= 0 || e < 0 || e >= 0.99
    bad = true;
end

end
