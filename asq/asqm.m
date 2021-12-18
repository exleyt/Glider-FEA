function [m,fm,whole] = asqm(R,sx3,K,a,fa,b,fb)
% Splits green integrand from a to b in half and calculates its asq value
% 
% Takes the green integrand, bounds a and b, and the integrand evaluated 
% at those bounds, fa and fb
% Returns the midpoint of a and b, m and the integrand evalulated at m, fm
% Also returns the asq integral value from a to b, whole
    m = (a + b)/2;
    fm = besselj(0,m*R)*exp(m*sx3)/(m-K);
    whole = abs(b - a)*(fa + 4*fm + fb)/6;
end

