function [result] = adaptive_simpson_quadrature_recursion(f,a,fa,b,fb,m,fm,whole,e)
%ADAPTIVE_SIMPSON_QUADRATURE_RECURSION Summary of this function goes here
%   Detailed explanation goes here
    [lm,flm,left] = adaptive_simpson_quadrature_mem(f,a,fa,m,fm);
    [rm,frm,right] = adaptive_simpson_quadrature_mem(f,m,fm,b,fb);
    delta = left + right - whole;
    if abs(delta) <= 15 * e
        result = left + right + delta / 15;
    else
        result = adaptive_simpson_quadrature_recursion(f,a,fa,m,fm,lm,flm,left,e/2);
        result = result + adaptive_simpson_quadrature_recursion(f,m,fm,b,fb,rm,frm,right,e/2); 
    end
end

