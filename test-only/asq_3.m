function [result] = asq_3(f,K,u,a,b,e,dmax)
% Numerically estimates the integral of f from a to b, error e depth dmax
% 
% Recursively estimates the integral of a triple parameter function f
% Integrates from a to b over the third argument where 
% K is the constant first argument
% u is the constant second argument
% a <= b
% e is a mimimum error value. Smaller errors take longer
% dmax is the maximum search depth or number of recursions
    fa = f(K,u,a);
    fb = f(K,u,b);
    [m,fm,whole] = asqm_3(f,K,u,a,fa,b,fb);
    result = asqr_3(f,K,u,a,fa,m,fm,b,fb,whole,e,0,dmax);
end