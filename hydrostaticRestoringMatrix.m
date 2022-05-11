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
% mesh is a triangulation object of the model
    FN = -faceNormal(mesh);
    [S,Sx,Sy,Sxy,Sxx,Syy] = waterplaneMoments(mesh.ConnectivityList,mesh.Points,FN);
    [V,CB] = volumeMoments(model.Mesh);
    [m,Cg] = massMoments(pm);

    pg = p*g;
    W = m*g; % weight
    Fb = pg*V; % buoyancy force
    Fd = Fb*CB - W*Cg;

    C33 = pg*S;
    C34 = pg*Sy;
    C35 = -pg*Sx;
    C44 = pg*Syy + Fd(3);
    C45 = -pg*Sxy;
    C46 = -Fb*CB(1) + W*Cg(3);
    C55 = pg*Sxx + Fd(3);
    C56 = -Fd(2);

    C = [
        0,  0,  0,      0,      0,          0;
        0,  0,  0,      0,      0,          0;
        0,  0,  C33,    C34,    C35,        0;
        0,  0,  C34,    C44,    C45,        C46;
        0,  0,  C35,    C45,    C55,        C56;
        0,  0,  0,      0,      0,          0;
    ];
end