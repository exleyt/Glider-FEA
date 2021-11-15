function [result] = adaptive_simpson_quadrature(f,a,b,e,dmax)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    fa = f(a);
    fb = f(b);
    [m,fm,whole] = adaptive_simpson_quadrature_mem(f, a, fa, b, fb);
    result = adaptive_simpson_quadrature_recursion(f, a, fa, m, fm, b, fb, whole, e, 0, dmax);
end

