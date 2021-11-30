function [result] = asqr(f,a,fa,m,fm,b,fb,whole,e,d,dmax)
% Recursivley estimates the integral of f from a to b, error e depth dmax
% 
% Recursively estimes the integral of a single parameter function f
% Integrates from a to b where
% a <= b
% m is the midpoint between and b
% fa, fb, and fm are the values of f evaluated at a, b, and m
% whole is the asq integral value from a to b
% e is a mimimum error value. Smaller errors take longer
% d is the current recursive depth
% dmax is the maximum search depth or number of recursions
    % Gets the asq integral value of the left side of the bounds
    [lm,flm,left] = asqm(f,a,fa,m,fm);
    % Gets the asq integral value of the right side of the bounds
    [rm,frm,right] = asqm(f,m,fm,b,fb);
    % Difference between left + right asq integral values and whole
    delta = left + right - whole;
    % Returns a result if delta is less than the epsilon value
    if abs(delta) <= 15 * e
        result = left + right + delta / 15;
    % Returns a result if max depth is surpassed
    elseif d >= dmax
        result = left + right + delta / 15;
    else
        % Recurse for the left side of the bounds
        result = asqr(f,a,fa,lm,flm,m,fm,left,e/2,d+1,dmax);
        % Recurse for the right side of the bounds
        result = result + asqr(f,m,fm,rm,frm,b,fb,right,e/2,d+1,dmax); 
    end
end