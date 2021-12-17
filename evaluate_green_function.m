function [result] = evaluate_green_function(x,xi,K)
% Numerically calculates the Green Function hdl.handle.net/20.500.12489/844
%
% Finds the complex valued influence of a pulsating source located at xi on
% the jth potential at x where: 
% K=w^2/g
% w is the waves angular frequency
% g is acceleration due to gravity
    dx2 = (x - xi).^2;
    r1 = sqrt(dx2(1) + dx2(2) + dx2(3));
    sx3 = x(3) + xi(3);
    r2 = sqrt(dx2(1) + dx2(2) + sx3^2);
    R = sqrt(dx2(1) + dx2(2));
    fun = @(k) besselj(0,k*R)*exp(k*sx3)/ (k-K);
    e = 10^-2; % an error value for asq estimation.
    e2 = 10^-3; % an epsilon value for isolation of the singularity.
    infi = 10; % a maximum depth for asq estimation.
    d = 20; % an approximate infinite upper bound for the integral.
    result = 1/r1 + 1/r2 + 2*K*(asq(fun,0,K-e2,e,d) + ...
        asq(fun,K+e2,infi,e,d)) - 2*pi*1i*K*exp(K*(x(3) + ...
        xi(3)))*besselj(0,K*R);
end