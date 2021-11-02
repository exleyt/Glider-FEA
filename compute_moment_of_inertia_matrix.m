function [moiMatrix] = compute_moment_of_inertia_matrix()
% A 3x3 describing the body's moment of inertia
%
% The moment of inertia matrix is defined as:
% [
%  Ixx    -Ixy     -Ixz
% -Ixy     Iyy     -Iyz
% -Ixz    -Iyz      Izz
% ]
% Where 
% Ixx = Integral(y^2+z^2)dm
% Iyy = Integral(x^2+z^2)dm
% Izz = Integral(x^2+y^2)dm
% Ixy = Integral(xy)dm
% Ixz = Integral(xz)dm
% Iyz = Integral(yz)dm
% m is mass
moiMatrix = [1,2,3;4,5,6;7,8,9];
end

