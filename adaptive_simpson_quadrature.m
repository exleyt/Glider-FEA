function [result] = adaptive_simpson_quadrature(f,a,b,e)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    syms k
    fa = f(a);
    fb = f(b);
    [m,fm,whole] = adaptive_simpson_quadrature_mem(f, a, fa, b, fb);
    result = adaptive_simpson_quadrature_recursion(f, a, fa, b, fb, m, fm, whole, e);
end

