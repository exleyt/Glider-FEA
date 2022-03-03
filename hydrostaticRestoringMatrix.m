function [C] = hydrostaticRestoringMatrix(pm,g,p,model,face)
% A 6x6 matrix describing the body's hydrostatic restoring force
%
% The matrix C is made up of the terms where
% pointMasses is the matrix of [mjs,pjs]
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
    [CL,P] = waterplaneTriangulation(model.Mesh,face);
    [S,CF,Ixx,Iyy] = waterplaneMoments(CL,P);
    [V,CB] = volumeMoments(model.Mesh);

    m = sum(pm(:,1)); % total mass
    Cg = sum(pm(:,2:4).*pm(:,1),1)/m; % center of gravity

    Ch = p*g*S; % force per unit distance
    W = m*g; % weight
    Fb = p*g*V; % buoyancy force
    C44 = -W*Cg(3) + Fb*CB(3) + p*g*Ixx + Ch*CF(2)^2;
    C55 = -W*Cg(3) + Fb*CB(3) + p*g*Iyy + Ch*CF(1)^2;

    C = [
        0,  0,  0,          0,                  0,                  0;
        0,  0,  0,          0,                  0,                  0;
        0,  0,  Ch,         Ch*CF(2),           -Ch*CF(1),          0;
        0,  0,  Ch*CF(2),   C44,                -Ch*CF(1)*CF(2),    (W*Cg(1) - Fb*CB(1));
        0,  0,  -Ch*CF(1),  -Ch*CF(1)*CF(2),    C55,                (W*Cg(2) - Fb*CB(2));
        0,  0,  0,          0,                  0,                  0;
    ];
end