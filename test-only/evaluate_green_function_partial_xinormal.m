function [df] = evaluate_green_function_partial_xinormal(x,xi,K,e,d,n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    addpath('..')
    grad = zeros(3,1);
    f = evaluate_green_function(x,xi,K,e,d);
    e2 = e^2;
    grad(1) = evaluate_green_function_partial_xi(x,xi,K,e2,d,1,f);
    grad(2) = evaluate_green_function_partial_xi(x,xi,K,e2,d,2,f);
    grad(3) = evaluate_green_function_partial_xi(x,xi,K,e2,d,3,f);
    df = dot(grad, n);
end

