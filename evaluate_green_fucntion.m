function [result] = evaluate_green_fucntion(x,y,z,rx,ry,rz,g,w)
%CALCULATE_GREEN_FUCNTION Summary of this function goes here
%   Detailed explanation goes here
    syms PV k
    K = w^2/g;
    r1 = sqrt((x - rx)^2 + (y - ry)^2 + (z - rz)^2);
    r2 = sqrt((x - rx)^2 + (y - ry)^2 + (z + rz)^2);
    R = sqrt((x - rx)^2 + (y - ry)^2);
    fun(k) = 1/(k-K)*exp(k*(z+rz))*besselj(0,k*R);
    maxk = log(realmax/(z+rz));
    result = 1/r1 + 1/r2 + 2*K*1*(adaptive_simpson_quadrature(fun,0,maxk,0.1)) - 2*pi*1i*K*exp(K*(z + rz))*besselj(0,K*R);
end