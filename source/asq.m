function [result] = asq(R,sx3,K,a,b)
% Numerically estimates the integral of the green integrand from a to b
% 
% Recursively estimates the integral of besselj(0,k*R)*exp(k*sx3)/(k-K)
% Integrates k from a to b where
% a <= b
    fa = besselj(0,a*R)*exp(a*sx3)/(a-K);
    fb = besselj(0,b*R)*exp(b*sx3)/(b-K);

    % Gets the asq integral value of the whole of the bounds
    m = (a + b)/2;
    fm = besselj(0,m*R)*exp(m*sx3)/(m-K);
    whole = (b - a)*(fa + 4*fm + fb)/6;
    
    e = 10^-4; % an error value for asq estimation.
    d = 30; % a maximum depth for asq estimation.
    result = asqr(R,sx3,K,a,b,fa,fm,fb,whole,e,d);
end