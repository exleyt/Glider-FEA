function [result] = evaluate_green_function_partial_xi(x,xi,K,j,f)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    e = 10^-2;
    xi(j) = xi(j) - e;
    fxix = evaluate_green_function(x,xi,K);
    result = (fxix - f) / e;
end

