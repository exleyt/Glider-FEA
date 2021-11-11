function [result] = adaptive_simpson_quadrature(f,a,b)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    fa = limit(f(k),k,a);
    fb = limit(f(k),k,b);
    [m,fm,whole] = adaptive_simpson_quadrature_mem(f, a, fa, b, fb);
    result = adaptive_simpson_quadrature_recursion(f, a, fa, b, fb, m, fm, whole);
end

