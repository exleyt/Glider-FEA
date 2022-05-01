function [result] = asqr(R,sx3,K,a,b,fa,fm,fb,whole,e,d)
% Recursivley estimates the green integrand from a to b
% 
% Recursively estimes the integral of besselj(0,k*R)*exp(k*sx3)/(k-K)
% Integrates k from a to b where
% a <= b
% m is the midpoint between and b
% fa, fb, and fm are the values of the green integrand at a, b, and m
% whole is the asq integral value from a to b
% e is a mimimum error value for the current depth
% d is 1 - the current recursive depth
    m = (a + b)/2;
    h = (b - a)/2;

    % Gets the asq integral value of the left side of the bounds
    lm = (a + m)/2;
    flm = besselj(0,lm*R)*exp(lm*sx3)/(lm-K);
    left = (h/6)*(fa + 4*flm + fm);

    % Gets the asq integral value of the right side of the bounds
    rm = (m + b)/2;
    frm = besselj(0,rm*R)*exp(rm*sx3)/(rm-K);
    right = (h/6)*(fm + 4*frm + fb);

    % Difference between left + right asq integral values and whole
    delta = left + right - whole;

    % Returns a result if delta is less than fifteen times the epsilon 
    %  value or max depth is surpassed
    if abs(delta) <= 15*e || d < 1
        result = left + right + delta / 15;
    else
        % Recurse for the left and then right sides of the bounds
        result = asqr(R,sx3,K,a,m,fa,flm,fm,left,e/2,d - 1) + ...
                 asqr(R,sx3,K,m,b,fm,frm,fb,right,e/2,d - 1); 
    end
end