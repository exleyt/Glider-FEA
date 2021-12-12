function [f,df] = evaluate_green_function_and_partial_xinormal(x,xi,n,K,e,d)
%EVALUATE_GREEN_FUNCTION_NORMAL_DERIVATIVE Summary of this function goes here
%   Detailed explanation goes here
    grad = zeros(3,1);
    f = evaluate_green_function(x,xi,K,e,d);
    for j=1:3
        grad(j) = evaluate_green_function_partial_xi(x,xi,K,e,d,j,f);
    end
    df = dot(grad, n);
end

