function [C,S] = stumpff(z)
%STUMPFF Compute Stumpff functions C(z), S(z).

tol = 1e-12;
if z > tol
    sz = sqrt(z);
    C = (1 - cos(sz)) / z;
    S = (sz - sin(sz)) / (sz^3);
elseif z < -tol
    sz = sqrt(-z);
    C = (cosh(sz) - 1) / (-z);
    S = (sinh(sz) - sz) / (sz^3);
else
    C = 1/2;
    S = 1/6;
end
end
