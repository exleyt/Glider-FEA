function [f,df] = evaluate_green_function_and_partial_xinormal(x,xi,n,K)
%EVALUATE_GREEN_FUNCTION_NORMAL_DERIVATIVE Summary of this function goes here
%   Detailed explanation goes here
    grad = zeros(3,1);
    f = evaluate_green_function(x,xi,K);
    for j=1:3
        grad(j) = evaluate_green_function_partial_xi(x,xi,K,j,f);
    end
    df = dot(grad, n);
end

