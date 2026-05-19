function state = KM_solver(t_sec, body)
% KM_solver - Keplerian body ephemeris in Mythara inertial frame.
% state = [r; v] in km and km/s.

mu = 139348062043.343;

a = body.a;
e = body.e;

if isfield(body,'i_rad')
    inc = body.i_rad;
else
    inc = deg2rad(body.i);
end

if isfield(body,'raan_rad')
    RAAN = body.raan_rad;
elseif isfield(body,'omega_rad')
    RAAN = body.omega_rad;
elseif isfield(body,'raan')
    RAAN = deg2rad(body.raan);
else
    RAAN = deg2rad(body.omega);
end

if isfield(body,'argp_rad')
    argp = body.argp_rad;
elseif isfield(body,'w_rad')
    argp = body.w_rad;
elseif isfield(body,'argp')
    argp = deg2rad(body.argp);
else
    argp = deg2rad(body.w);
end

if isfield(body,'M0_rad')
    M0 = body.M0_rad;
else
    M0 = deg2rad(body.M0);
end

n = sqrt(mu/a^3);
M = M0 + n*t_sec;
M = mod(M,2*pi);

E = solve_kepler(M,e);

cE = cos(E);
sE = sin(E);

r_pf = [a*(cE - e);
        a*sqrt(1-e^2)*sE;
        0];

r = a*(1 - e*cE);

v_pf = sqrt(mu*a)/r * [-sE;
                        sqrt(1-e^2)*cE;
                        0];

R3O = [ cos(RAAN) -sin(RAAN) 0;
        sin(RAAN)  cos(RAAN) 0;
        0          0         1];

R1i = [1 0 0;
       0 cos(inc) -sin(inc);
       0 sin(inc)  cos(inc)];

R3w = [ cos(argp) -sin(argp) 0;
        sin(argp)  cos(argp) 0;
        0          0         1];

Q = R3O*R1i*R3w;

state = [Q*r_pf; Q*v_pf];

end

function E = solve_kepler(M,e)

E = M;
if e > 0.8
    E = pi;
end

for k = 1:100
    dE = -(E - e*sin(E) - M)/(1 - e*cos(E));
    E = E + dE;
    if abs(dE) < 1e-14
        return
    end
end

error('Kepler solver did not converge')

end