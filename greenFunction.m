function [result] = greenFunction(x,xi,K)
% Numerically calculates the Green Function for finding the RAO
%
% Finds the complex valued influence of a pulsating source located at xi on
% the jth potential at x where: 
% x is the location of the potential
% xi is the location of the pulsating source
% K = w^2/g
% w is the waves angular frequency
% g is acceleration due to gravity
    dx2 = (x - xi).^2;
    r1 = sqrt(dx2(1) + dx2(2) + dx2(3));
    sx3 = x(3) + xi(3);
    r2 = sqrt(dx2(1) + dx2(2) + sx3^2);
    R = sqrt(dx2(1) + dx2(2));
    e = 10^-3; % an epsilon value for isolation of the singularity.
    infi = 10; % an approximate infinite upper bound for the integral.
    result = 1/r1 + 1/r2 + 2*K*(asq(R,sx3,K,0,K-e) + ...
        asq(R,sx3,K,K+e,infi)) - 2*pi*1i*K*exp(K*(x(3) + ...
        xi(3)))*besselj(0,K*R);
end