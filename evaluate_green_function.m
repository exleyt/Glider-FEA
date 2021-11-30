function [result] = evaluate_green_function(x,xi,K)
% Numerically calculates the Green Function hdl.handle.net/20.500.12489/844
%
% Finds the complex valued influence of a pulsating source located at xi on
% the jth potential at x where: 
% K=w^2/g
% w is the waves angular frequency
% g is acceleration due to gravity
    % An error value for asq. Needs optimization! Too low and all values
    % become 0. Too high and results will suffer.
    e = 10 ^ -6; 
    % An asq estimation of infinity. Should be around the value that fun->0
    infi = 20; 
    % A maximum depth for asq estimation.
    d = 20;
    r1 = sqrt((x(1) - xi(1))^2 + (x(2) - xi(2))^2 + (x(3) - xi(3))^2);
    r2 = sqrt((x(1) - xi(1))^2 + (x(2) - xi(2))^2 + (x(3) + xi(3))^2);
    R = sqrt((x(1) - xi(1))^2 + (x(2) - xi(2))^2);
    fun = @(k) besselj(0,k*R)*exp(k*(x(3)+xi(3)))/(k-K);
    result = 1/r1 + 1/r2 + 2*K*1*(asq(fun,0,K-e,e,d) + asq(fun,K+e,K+e+1,e,d) + asq(fun,K+e+1,infi,e,d)) - 2*pi*1i*K*exp(K*(x(3) + xi(3)))*besselj(0,K*R);
end