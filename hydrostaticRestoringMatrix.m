function [C] = hydrostaticRestoringMatrix(pm,g,p,model,mesh)
% A 6x6 matrix describing the body's hydrostatic restoring force
%
% The matrix C is made up of the terms where
% pm is the matrix of [mjs,pjs]
% mj is the jth point mass
% pj is the position vector of the jth point mass (xj,yj,zj)
% g is the acceleration due to gravity
% p is the water density
% model is a meshed pdemodel
% face is the face number of the waterplane on the model's geometry
% center of gravity, Cg is [xg,yg,zg]
% center of buoyancy, Cb is [xb,yb,zb]
% m is the total glider's mass
%
% This is from a textbook where it is assumed that the origin is the center
%  of flotation which it very much is not as currently defined
    %[CL,P] = waterplaneTriangulation(model.Mesh,face);
    FN = -faceNormal(mesh);
    [S,Sx,Sy,Sxy,Sxx,Syy] = waterplaneMoments(mesh.ConnectivityList,mesh.Points,FN);
    [V,CB] = volumeMoments(model.Mesh);
    [m,Cg] = massMoments(pm);

    pg = p*g;
    W = m*g; % weight
    Fb = pg*V; % buoyancy force
    Fd = Fb*CB - W*Cg;
    C44 = pg*Syy + Fd(3);
    C55 = pg*Sxx + Fd(3);

    C = [
        0,  0,  0,      0,      0,          0;
        0,  0,  0,      0,      0,          0;
        0,  0,  pg*S,   pg*Sy,  -pg*Sx,     0;
        0,  0,  pg*Sy,      C44,    -pg*Sxy,    -Fb*CB(1) + W*Cg(3);
        0,  0,  -pg*Sx,      -pg*Sxy,      C55,        -Fd(2);
        0,  0,  0,      0,      0,          0;
    ];
end