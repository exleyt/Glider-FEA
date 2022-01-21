function [moiMatrix] = MomentOfInertia(pointMasses)
% A 3x3 describing the body's moment of inertia
%
% The moment of inertia matrix is the sum of the matrices: 
% [
%   Ixx     Ixy     Ixz
%   Ixy     Iyy     Iyz
%   Ixz     Iyz     Izz
% ]

% Where
% Ixx = Sum(1,j,(yj^2 + zj^2)mj)
% Iyy = Sum(1,j,(xj^2 + zj^2)mj)
% Izz = Sum(1,j,(xj^2 + yj^2)mj)
% Ixy = -Sum(1,j,xj * yj * mj)
% Ixz = -Sum(1,j,xj * zj * mj)
% Iyz = -Sum(1,j,yj * zj * mj)   
% pointMasses is the matrix of mjs and pjs
% mj is the jth point mass
% pj is the position vector of the jth point mass (xj,yj,zj)
%
% It is assumed that the glider can be treated as a sum of point masses and
% that pointMasses are have positions relative to the gliders center of
% mass
    xx = 0;
    yy = 0;
    zz = 0;
    xy = 0;
    xz = 0;
    yz = 0;
    for pm = pointMasses
        xx = xx + pm(1)*(pm(3)^2 + pm(4)^2);
        yy = yy + pm(1)*(pm(2)^2 + pm(4)^2);
        zz = zz + pm(1)*(pm(2)^2 + pm(3)^2);
        xy = xy + pm(1)*pm(2)*pm(3);
        xz = xz + pm(1)*pm(2)*pm(4);
        yz = yz + pm(1)*pm(3)*pm(4);
    end
    moiMatrix = [
        xx, -xy, -xz;
       -xy,  yy, -yz;
       -xz, -yz,  zz;
    ];
end
% TODO 
% Since the glider actually rotates about its center of buoyancy, pjs
% should probably be given in relation to that instead.
