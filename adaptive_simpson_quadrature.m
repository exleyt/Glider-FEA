function [result] = adaptive_simpson_quadrature(f,a,b,e,dmax)
% Numerically estimates the integral of f from a to b, error e depth dmax
% 
% Recursively estimes the integral of a single parameter function f
% Integrates from a to b where
% a <= b
% e is a mimimum error value. Smaller errors take longer
% dmax is the maximum search depth or number of recursions
    fa = f(a);
    fb = f(b);
    [m,fm,whole] = adaptive_simpson_quadrature_mem(f, a, fa, b, fb);
    result = adaptive_simpson_quadrature_recursion(f, a, fa, m, fm, b, fb, whole, e, 0, dmax);
end

