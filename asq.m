function [result] = asq(R,sx3,K,a,b)
% Numerically estimates the integral of the green integrand from a to b
% 
% Recursively estimates the integral of besselj(0,k*R)*exp(k*sx3)/(k-K)
% Integrates k from a to b where
% a <= b
    fa = besselj(0,a*R)*exp(a*sx3)/(a-K);
    fb = besselj(0,b*R)*exp(b*sx3)/(b-K);
    [m,fm,whole] = asqm(R,sx3,K,a,fa,b,fb);
    e = 10^-2; % an error value for asq estimation.
    result = asqr(R,sx3,K,a,fa,m,fm,b,fb,whole,e,0);
end