function [m,fm,whole] = asqm(f,a,fa,b,fb)
% Splits f from a to b down the middle and calculates its asq value
% 
% Takes the integrand function f, bounds a and b, and the function 
% evaluated at those bounds, fa and fb
% Returns the midpoint of a and b, m and the fucntion evalulated at m, fm
% Also returns the asq integral value from a to b, whole
    m = (a + b)/2;
    fm = f(m);
    whole = abs(b - a)*(fa + 4*fm + fb)/6;
end

