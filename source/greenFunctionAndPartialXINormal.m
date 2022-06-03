function [f,df] = greenFunctionAndPartialXINormal(x,xi,FN,K)
% Returns the green function and its partial normal of xi
%
% Finds the complex valued influence of a pulsating source located at xi on
%  the potential at x and its normal derivative where:
% x is the location of the potential
% xi is the location of the pulsating source
% FN is the normal vector of the tile corresponding to xi 
% K = w^2/g
    f = greenFunction(x,xi,K); % the green function evaluated at xi
    e = 10^-8; % epsilon for estimating 
    fe = greenFunction(x,xi + e*FN,K);
    df = (fe - f)/e; 
end