function [result] = greenFunctionPartialXI(x,xi,K,j,f)
% Returns the green function and its partial normal of xi
%
% Finds the partial derivative of the complex valued influence of a 
%  pulsating source located at xi w.r.t. xi(j) where:
% x is the location of the potential
% xi is the location of the pulsating source
% K = w^2/g
% w is the waves angular frequency
% g is acceleration due to gravity
    e = 10^-2; % epsilon for estimating 
    xi(j) = xi(j) + e;
    fe = greenFunction(x,xi,K);
    result = (fe - f) / e;
end

