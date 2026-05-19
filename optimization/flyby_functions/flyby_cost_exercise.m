function [cost, met] = flyby_cost_exercise(x, v_in, t_fb, id_fb, id_next, bodies)
% Cost for closing one patched-conics flyby followed by one Lambert leg.
% x = [TOF_years, Lambert_path]

met = struct();
met.valid = false;
met.mag_err = NaN;
met.h = NaN;
met.h_over_R = NaN;
met.delta_deg = NaN;
met.vinf_in_mag = NaN;
met.vinf_out_mag = NaN;
met.vinf_arr_mag = NaN;

cost = 1e6;

y2sec = bodies.year_sec;
mu_star = bodies.mu_star;

TOF = x(1);
path = round(x(2));

if ~isfinite(TOF) || TOF <= 0 || ~ismember(path, [0 1])
    return;
end

try
    tA = t_fb*y2sec;
    tB = (t_fb + TOF)*y2sec;

    state_fb   = KM_solver(tA, bodies.(id_fb));
    state_next = KM_solver(tB, bodies.(id_next));

    state_fb   = state_fb(:);
    state_next = state_next(:);

    r_fb   = state_fb(1:3);
    v_p_fb = state_fb(4:6);
    r_next = state_next(1:3);
    v_next = state_next(4:6);

    % Correct signature:
    % [v1,v2,iter] = lambert_universal(mu, t1, r1, t2, r2, path)
    [v_dep, v_arr] = lambert_universal(mu_star, tA, r_fb, tB, r_next, path);

    v_dep = v_dep(:);
    v_arr = v_arr(:);

    v_inf_in  = v_in(:) - v_p_fb;
    v_inf_out = v_dep   - v_p_fb;
    v_inf_arr = v_arr   - v_next;

    vin_mag_in  = norm(v_inf_in);
    vin_mag_out = norm(v_inf_out);

    met.vinf_in_mag  = vin_mag_in;
    met.vinf_out_mag = vin_mag_out;
    met.vinf_arr_mag = norm(v_inf_arr);
    met.mag_err = abs(vin_mag_out - vin_mag_in);

    if vin_mag_in <= 0 || vin_mag_out <= 0
        return;
    end

    cos_delta = dot(v_inf_in, v_inf_out)/(vin_mag_in*vin_mag_out);
    cos_delta = max(-1, min(1, cos_delta));
    delta = acos(cos_delta);
    met.delta_deg = rad2deg(delta);

    body = bodies.(id_fb);

    % Nyxar / massless case
    if body.mu <= 0
        dir_err = norm(v_inf_out - v_inf_in);
        cost = 1e8*met.mag_err^2 + 1e8*dir_err^2;
        met.valid = (met.mag_err < 1e-7) && (dir_err < 1e-7);
        return;
    end

    if delta <= 0
        return;
    end

    rp = body.mu/vin_mag_in^2 * (1/sin(delta/2) - 1);
    h  = rp - body.radius;

    met.h = h;
    met.h_over_R = h/body.radius;

    alt_penalty = 0;
    if met.h_over_R < 0.1
        alt_penalty = alt_penalty + (0.1 - met.h_over_R)^2;
    elseif met.h_over_R > 100
        alt_penalty = alt_penalty + (met.h_over_R - 100)^2;
    end

    % Main objective: close v_inf magnitude; keep altitude feasible.
    cost = 1e8*met.mag_err^2 + 1e4*alt_penalty + 1e-3*met.vinf_arr_mag^2;

    met.valid = (met.mag_err < 1e-7) && ...
                (met.h_over_R >= 0.1) && ...
                (met.h_over_R <= 100);

catch
    cost = 1e6;
end

end
