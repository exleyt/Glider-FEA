function [f,df] = greenFunctionAndPartialXINormal(x,xi,n,K)
%EVALUATE_GREEN_FUNCTION_NORMAL_DERIVATIVE Summary of this function goes here
%   Detailed explanation goes here
    grad = zeros(3,1);
    f = greenFunction(x,xi,K);
    for j=1:3
        grad(j) = greenFunctionPartialXI(x,xi,K,j,f);
    end
    df = dot(grad, n);
end

