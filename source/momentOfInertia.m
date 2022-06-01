function [result] = momentOfInertia(pm)
% A 3x3 describing the body's moment of inertia
%
% The moment of inertia matrix is the sum of the matrices: 
% [
%   Ixx     Ixy     Ixz
%   Ixy     Iyy     Iyz
%   Ixz     Iyz     Izz
% ] s.t.
% Ixx = Sum(mj*(yj^2 + zj^2))
% Iyy = Sum(mj*(xj^2 + zj^2))
% Izz = Sum(mj*(xj^2 + yj^2))
% Ixy = -Sum(mj*xj*yj)
% Ixz = -Sum(mj*xj*zj)
% Iyz = -Sum(mj*yj*zj)   
% pm is the matrix of [mjs,pjs]
% mj is the jth point mass
% pj is the position vector of the jth point mass (xj,yj,zj)
% It is assumed that the glider can be treated as a sum of point masses
    Ixx = 0;
    Iyy = 0;
    Izz = 0;
    Ixy = 0;
    Ixz = 0;
    Iyz = 0;
    [N,~] = size(pm);
    for j = 1:N
        Ixx = Ixx + pm(j,1)*(pm(j,3)^2 + pm(j,4)^2);
        Iyy = Iyy + pm(j,1)*(pm(j,2)^2 + pm(j,4)^2);
        Izz = Izz + pm(j,1)*(pm(j,2)^2 + pm(j,3)^2);
        Ixy = Ixy + pm(j,1)*pm(j,2)*pm(j,3);
        Ixz = Ixz + pm(j,1)*pm(j,2)*pm(j,4);
        Iyz = Iyz + pm(j,1)*pm(j,3)*pm(j,4);
    end
    result = [
        Ixx, -Ixy, -Ixz;
       -Ixy,  Iyy, -Iyz;
       -Ixz, -Iyz,  Izz;
    ];
end