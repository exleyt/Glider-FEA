function [C] = hydrostaticRestoringMatrix(pm,g,p,mesh,face)
% A 6x6 matrix describing the body's hydrostatic restoring force
%
% The matrix C is made up of the terms, Cij s.t.
% C22 = p*g*S
% C44 = g*(p*(S11 + V*zb) - m*zg)
% C45 = -g*(p*V*yb - m*yg)
% C65 = -g*(p*V*xb - m*xg)
% C66 = g*(p*(S22 + V*zb) - m*zg)
% p is the water density
% g is the acceleration due to gravity
% V is the displaced volume
% center of gravity, cog is [xg,yg,zg]
% center of buoyancy, cob is [xb,yb,zb]
% m is the total glider's mass
%
% This is from a textbook where it is assumed that the origin is the center
%  of flotation which it very much is not as currently defined
    [CL,P] = waterplaneTriangulation(mesh,face);
    [S,CF,Ixx,Iyy] = waterplaneMoments(CL,P);
    [V,CB] = volumeMoments(app.Model.Mesh);

    [N,~] = size(pm);

    % Center of gravity
    CG = sum(pm(1:N,2:4),1) / N;

    % Mass of glider
    m = sum(pm(1:N,1));

    Ch = p*g*S;
    W = m*g;
    Fb = p*g*V;
    C44 = -W*CG(3) + Fb*CB(3) + p*g*Ixx + Ch*CF(2)^2;
    C55 = -W*CG(3) + Fb*CB(3) + p*g*Iyy + Ch*CF(1)^2;

    C = [
        0,  0,  0,          0,                  0,                  0;
        0,  0,  0,          0,                  0,                  0;
        0,  0,  Ch,         Ch*CF(2),           -Ch*CF(1),          0;
        0,  0,  Ch*CF(2),   C44,                -Ch*CF(1)*CF(2),    (W - Fb)*CF(1);
        0,  0,  -Ch*CF(1),  -Ch*CF(1)*CF(2),    C55,                (W - Fb)*CF(2);
        0,  0,  0,          0,                  0,                  0;
    ];
end

