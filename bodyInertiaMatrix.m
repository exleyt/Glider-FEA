function [bicMatrix] = bodyInertiaMatrix(pm)
% A 6x6 matrix describing the body's inertia  with 6 degrees of freedom
% 
% The body's inertia coefficients matrix is defined as:
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
    [~,N] = size(pm);     
    cog = zeros(3,1); 
    for j = 1:3
        cog(j) = sum(pm(j+1,1:N)) / N;
    end
    m = sum(pm(1,1:N));
    M = diag([m,m,m]);
    A = [
        0, m*cog(3), -m*cog(2);
        -m*cog(3), 0, m*cog(1);
        m*cog(2), -m*cog(1), 0;
    ];
    I = MomentOfInertia(pm);
    bicMatrix = [  
        M, A;
        -A, I;
    ];
end
% TODO
% Make compute_moment_of_inertia_matrix()
% Update compute_moment_of_inertia_matrix arguments
% Update compute_body_inertia_coefficients_matrix arguments
