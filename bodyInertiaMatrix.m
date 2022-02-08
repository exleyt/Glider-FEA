function [result] = bodyInertiaMatrix(pm)
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
% ] Where:
% pointMasses is the matrix of [mjs,pjs]
% mj is the jth point mass
% pj is the position vector of the jth point mass (xj,yj,zj)
% the center of gravity vector (cog) is [xg, yg, zg]
% I is the 3x3 moment of inertia matrix
    [j,~] = size(pm);     
    
    cog = zeros(3,1); 
    for i = 1:3
        cog(i) = sum(pm(1:j,i+1)) / j;
    end
    
    m = sum(pm(1:j,1));

    M = diag([m,m,m]);

    A = [
        0, m*cog(3), -m*cog(2);
        -m*cog(3), 0, m*cog(1);
        m*cog(2), -m*cog(1), 0;
    ];

    I = momentOfInertia(pm);

    result = [  
        M, A;
        -A, I;
    ];
end