function [result] = evaluate_green_fucntion(x,y,z,xix,xiy,xiz,K)
% Numerically calculates the Green Function hdl.handle.net/20.500.12489/844
%
% Finds the complex valued influence of a pulsating source located at xi on
% the jth potential at x where: 
% K=w^2/g
% w is the waves angular frequency
% g is acceleration due to gravity
    % An error value for asq. Needs optimization! Too low and all values
    % become 0. Too high and results will suffer.
    e = 10 ^ -16; 
    % An asq estimation of infinity. Should be around the value that fun->0
    infi = 100; 
    r1 = sqrt((x - xix)^2 + (y - xiy)^2 + (z - xiz)^2);
    r2 = sqrt((x - xix)^2 + (y - xiy)^2 + (z + xiz)^2);
    R = sqrt((x - xix)^2 + (y - xiy)^2);
    fun = @(k) besselj(0,k*R)*exp(k*(z+xiz))/(k-K);
    result = 1/r1 + 1/r2 + 2*K*1*(adaptive_simpson_quadrature(fun,0,K-e,10^-6,30) + adaptive_simpson_quadrature(fun,K+e,K+e+1,10^-6,30) + adaptive_simpson_quadrature(fun,K+e+1,infi,10^-6,30)) - 2*pi*1i*K*exp(K*(z + xiz))*besselj(0,K*R);
end