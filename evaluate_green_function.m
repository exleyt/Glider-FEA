function [result] = evaluate_green_function(x,xi,K,e,d)
% Numerically calculates the Green Function hdl.handle.net/20.500.12489/844
%
% Finds the complex valued influence of a pulsating source located at xi on
% the jth potential at x where: 
% K=w^2/g
% w is the waves angular frequency
% g is acceleration due to gravity
% e is an error value for asq estimation. Needs optimization!
% d is a maximum depth for asq estimation. Needs optimization!
    % An asq estimation of infinity. Should be around the value that fun->0
    % Needs optimization, but not essential. It saves very little time.
    infi = 100;
    dx2 = (x - xi)^2;
    r1 = sqrt(dx2(1) + dx2(2) + dx2(3));
    r2 = sqrt(dx2(1) + dx2(2) + (x(3) + xi(3))^2);
    R = sqrt(dx2(1) + dx2(2));
    fun = @(k) besselj(0,k*R)*exp(k*(x(3)+xi(3)))/(k-K);
    result = 1/r1 + 1/r2 + 2*K*1*(asq(fun,0,K-e,e,d) + asq(fun,K+e,infi,e,d)) - 2*pi*1i*K*exp(K*(x(3) + xi(3)))*besselj(0,K*R);
end