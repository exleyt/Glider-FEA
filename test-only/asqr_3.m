function [result] = asqr_3(f,K,u,a,fa,m,fm,b,fb,whole,e,d,dmax)
% Recursivley estimates the integral of f from a to b, error e depth dmax
% 
% Recursively estimes the integral of a triple parameter function f
% Integrates from a to b where over the third argument where
% K is the constant first argument
% u is the constant second argument
% a <= b
% m is the midpoint between and b
% fKa, fKb, and fKm are the values of f evaluated at a, b, and m
% whole is the asq integral value from a to b
% e is a mimimum error value. Smaller errors take longer
% d is the current recursive depth
% dmax is the maximum search depth or number of recursions
    % Gets the asq integral value of the left side of the bounds
    [lm,flm,left] = asqm_3(f,K,u,a,fa,m,fm);
    % Gets the asq integral value of the right side of the bounds
    [rm,frm,right] = asqm_3(f,K,u,m,fm,b,fb);
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
        result = asqr_3(f,K,u,a,fa,lm,flm,m,fm,left,e/2,d+1,dmax);
        % Recurse for the right side of the bounds
        result = result + asqr_3(f,K,u,m,fm,rm,frm,b,fb,right,e/2,d+1,dmax); 
    end
end