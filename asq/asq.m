function [result] = asq(f,a,b,e,dmax)
% Numerically estimates the integral of f from a to b, error e depth dmax
% 
% Recursively estimates the integral of a single parameter function f
% Integrates from a to b where
% a <= b
% e is a mimimum error value. Smaller errors take longer
% dmax is the maximum search depth or number of recursions
    fa = f(a);
    fb = f(b);
    [m,fm,whole] = asqm(f,a,fa,b,fb);
    result = asqr(f,a,fa,m,fm,b,fb,whole,e,0,dmax);
end