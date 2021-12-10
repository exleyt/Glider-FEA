function [result] = evaluate_green_function_partial_xi(x,xi,K,e,d,j,f)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    xi(j) = xi(j) - e;
    fxix = evaluate_green_function(x,xi,K,e,d);
    result = (fxix - f) / e;
end

