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

    C = zeros(6,6);
    C(3,3) = p*g*S;
    C(4,4) = g*(p*(Iyy + V*CB(3)) - m*CG(3));
    C(4,6) = -g*(p*V*CB(1) - m*CG(1));
    C(5,5) = g*(p*(Ixx + V*CB(3)) - m*CG(3));
    C(5,6) = -g*(p*V*CB(2) - m*CG(2));
end

