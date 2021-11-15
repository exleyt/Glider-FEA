function [m,fm,whole] = adaptive_simpson_quadrature_mem(f,a,fa,b,fb)
%ADAPTIVE_SIMPSON_QUADRATURE_MEM Summary of this function goes here
%   Detailed explanation goes here
    m = (a + b)/2;
    fm = f(m);
    whole = abs(b - a)*(fa + 4*fm + fb)/6;
end

