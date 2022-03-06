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
% pm is the matrix of [mjs,pjs]
% mj is the jth point mass
% pj is the position vector of the jth point mass (xj,yj,zj)
% the center of gravity vector (Cg) is [xg, yg, zg]
% I is the 3x3 moment of inertia matrix
    [m,Cg] = massMoments(pm);
    
    M = diag([m,m,m]);

    A = [
        0, m*Cg(3), -m*Cg(2);
        -m*Cg(3), 0, m*Cg(1);
        m*Cg(2), -m*Cg(1), 0;
    ];

    I = momentOfInertia(pm);

    result = [  
        M, A;
        -A, I;
    ];
end