function [result] = evaluate_green_fucntion(x,y,z,rx,ry,rz,K)
%CALCULATE_GREEN_FUCNTION Summary of this function goes here
%   Detailed explanation goes here
    e = 10 ^ -16;
    infi = 100;
    r1 = sqrt((x - rx)^2 + (y - ry)^2 + (z - rz)^2);
    r2 = sqrt((x - rx)^2 + (y - ry)^2 + (z + rz)^2);
    R = sqrt((x - rx)^2 + (y - ry)^2);
    fun = @(k) besselj(0,k*R)*exp(k*(z+rz))/(k-K);
    result = 1/r1 + 1/r2 + 2*K*1*(adaptive_simpson_quadrature(fun,0,K-e,10^-6,30) + adaptive_simpson_quadrature(fun,K+e,K+e+1,10^-6,30) + adaptive_simpson_quadrature(fun,K+e+1,infi,10^-6,30)) - 2*pi*1i*K*exp(K*(z + rz))*besselj(0,K*R);
end