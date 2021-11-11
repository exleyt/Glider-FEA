function [result] = evaluate_green_fucntion(x,y,z,rx,ry,rz,g,w)
%CALCULATE_GREEN_FUCNTION Summary of this function goes here
%   Detailed explanation goes here
    K = w^2/g;
    r1 = sqrt((x - rx)^2 + (y - ry)^2 + (z - rz)^2);
    r2 = sqrt((x - rx)^2 + (y - ry)^2 + (z + rz)^2);
    R = sqrt((x - rx)^2 + (y - ry)^2);
    result = 1/r1 + 1/r2 - 2*pi*1i*K*exp(K*(z + rz))*besselj(0,K*R);
end

