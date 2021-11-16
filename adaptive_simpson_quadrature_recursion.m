function [result] = adaptive_simpson_quadrature_recursion(f,a,fa,m,fm,b,fb,whole,e,d,dmax)
%ADAPTIVE_SIMPSON_QUADRATURE_RECURSION Summary of this function goes here
%   Detailed explanation goes here
    [lm,flm,left] = adaptive_simpson_quadrature_mem(f,a,fa,m,fm);
    [rm,frm,right] = adaptive_simpson_quadrature_mem(f,m,fm,b,fb);
    delta = left + right - whole;
    if abs(delta) <= 15 * e
        result = left + right + delta / 15;
    elseif d >= dmax
        result = left + right + delta / 15;
    else
        result = adaptive_simpson_quadrature_recursion(f,a,fa,lm,flm,m,fm,left,e/2,d+1,dmax);
        result = result + adaptive_simpson_quadrature_recursion(f,m,fm,rm,frm,b,fb,right,e/2,d+1,dmax); 
    end
end

