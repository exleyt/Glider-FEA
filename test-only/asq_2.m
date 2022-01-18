function [result] = asq_2(f,K,a,b,e,dmax)
% Numerically estimates the integral of f from a to b, error e depth dmax
% 
% Recursively estimates the integral of a double parameter function f
% Integrates from a to b over the second argument where 
% K is the constant first argument
% a <= b
% e is a mimimum error value. Smaller errors take longer
% dmax is the maximum search depth or number of recursions
    fa = f(K,a);
    fb = f(K,b);
    [m,fKm,whole] = asqm_2(f,K,a,fa,b,fb);
    result = asqr_2(f,K,a,fa,m,fKm,b,fb,whole,e,0,dmax);
end