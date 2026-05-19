function [v1,v2,iter] = lambert_universal(mu,t1,r1,t2,r2,path,Nz)
%LAMBERT_UNIVERSAL Robust universal-variable Lambert solver.
%
% Backward compatible:
%   [v1,v2,iter] = lambert_universal(mu,t1,r1,t2,r2,path)
% uses Nz = 3000, same robust bracket density as the original version.
%
% Fast search mode:
%   [v1,v2,iter] = lambert_universal(mu,t1,r1,t2,r2,path,Nz)
% use e.g. Nz = 400 or 500 inside the beam search.

if nargin < 7 || isempty(Nz)
    Nz = 3000;
end
Nz = max(50, round(Nz));

r1 = r1(:);
r2 = r2(:);

r1n = norm(r1);
r2n = norm(r2);
dt  = t2 - t1;

if dt <= 0
    error('Lambert: non-positive TOF')
end

cosd = dot(r1,r2)/(r1n*r2n);
cosd = max(-1,min(1,cosd));

dtheta = acos(cosd);
sind = sin(dtheta);

if path == 1
    sind = -sind;
end

A = sind*sqrt(r1n*r2n/(1 - cosd));

if abs(A) < 1e-12
    error('Lambert: singular geometry')
end

target = dt;
zmin = -4*pi^2 + 1e-8;
zmax =  80;

% Bracket scan. Nz is the main runtime knob.
zgrid = linspace(zmin,zmax,Nz);

zlo = NaN;
zhi = NaN;
Flo = NaN;

for k = 1:Nz
    z = zgrid(k);
    [F,ok] = time_eq(z,A,r1n,r2n,mu,target);

    if ~ok || ~isfinite(F)
        continue
    end

    if isnan(zlo)
        zlo = z;
        Flo = F;
        continue
    end

    if Flo*F <= 0
        zhi = z;
        break
    end

    zlo = z;
    Flo = F;
end

if isnan(zhi)
    error('Lambert: no z bracket found')
end

% Bisection. Keep Flo cached; do not recompute time_eq(zlo) every iteration.
for iter = 1:200
    z = 0.5*(zlo+zhi);
    [F,ok] = time_eq(z,A,r1n,r2n,mu,target);

    if ~ok || ~isfinite(F)
        zlo = z;
        Flo = F;
        continue
    end

    if abs(F) < 1e-10
        break
    end

    if Flo*F <= 0
        zhi = z;
    else
        zlo = z;
        Flo = F;
    end
end

[C,S] = stumpff(z);
y = r1n + r2n + A*(z*S - 1)/sqrt(C);

if y <= 0
    error('Lambert: negative y')
end

f  = 1 - y/r1n;
g  = A*sqrt(y/mu);
gd = 1 - y/r2n;

if abs(g) < 1e-12
    error('Lambert: singular g')
end

v1 = (r2 - f*r1)/g;
v2 = (gd*r2 - r1)/g;

end

function [F,ok] = time_eq(z,A,r1n,r2n,mu,target)

ok = false;

[C,S] = stumpff(z);

if C <= 0
    F = Inf;
    return
end

y = r1n + r2n + A*(z*S - 1)/sqrt(C);

if y <= 0
    F = Inf;
    return
end

x = sqrt(y/C);
tof = (x^3*S + A*sqrt(y))/sqrt(mu);

F = tof - target;
ok = true;

end
