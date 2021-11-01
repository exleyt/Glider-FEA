function [bicMatrix] = compute_body_inertia_coefficients_matrix(m,cog)
%M  A 6x6 matrix describing the body's inertia  with 6 degrees of freedom
% 
% The matrix of body inertia coefficients is defined as
% [
%   m       0       0       0       mzg     -myg
%   0       m       0       -mzg    0       mxg
%   0       0       m       myg     -mxg    0
%   0       -mzg    myg     I11     I12     I13
%   mzg     0       -mxg    I21     I22     I23
%   -myg    mgx     0       I31     I32     I33
% ]
% Where:
% m is the total body mass
% the center of gravity vector (cog) is [xg, yg, zg]
% I is the 3x3 moment of inertia matrix
M = diag([m,m,m]);
A = [
    0, m*cog(3), -m*cog(2);
    -m*cog(3), 0, m*cog(1);
    m*cog(2), -m*cog(1), 0;
];
I = compute_moment_of_inertia_matrix();
bicMatrix = [  
    M, A;
    -A, I;
];
end

