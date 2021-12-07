function [result] = evaluate_green_function_partial_xi(x,xi,K,e,d,j,e2,f)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    xi(j) = xi(j) - e2;
    fxix = evaluate_green_function(x,xi,K,e,d);
    result = (fxix - f) / e2;
end

