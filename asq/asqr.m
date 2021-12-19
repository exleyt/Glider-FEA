function [result] = asqr(R,sx3,K,a,fa,m,fm,b,fb,whole,e,d)
% Recursivley estimates the green integrand from a to b
% 
% Recursively estimes the integral of besselj(0,k*R)*exp(k*sx3)/(k-K)
% Integrates k from a to b where
% a <= b
% m is the midpoint between and b
% fa, fb, and fm are the values of the green integrand at a, b, and m
% whole is the asq integral value from a to b
% e is a mimimum error value for the current depth
% d is the current recursive depth
    % Gets the asq integral value of the left side of the bounds
    [lm,flm,left] = asqm(R,sx3,K,a,fa,m,fm);
    % Gets the asq integral value of the right side of the bounds
    [rm,frm,right] = asqm(R,sx3,K,m,fm,b,fb);
    % Difference between left + right asq integral values and whole
    delta = left + right - whole;
    % Returns a result if delta is less than the epsilon value
    if abs(delta) <= 15 * e
        result = left + right + delta / 15;
    % Returns a result if max depth is surpassed
    elseif d >= 20 % A maximum depth for asq estimation.
        result = left + right + delta / 15;
    else
        % Recurse for the left side of the bounds
        result = asqr(R,sx3,K,a,fa,lm,flm,m,fm,left,e/2,d+1);
        % Recurse for the right side of the bounds
        result = result + asqr(R,sx3,K,m,fm,rm,frm,b,fb,right,e/2,d+1); 
    end
end