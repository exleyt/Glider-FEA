function [C] = hydrostaticRestoringMatrix(pm,g,p,model)
% A 6x6 matrix describing the body's hydrostatic restoring force
%
% The matrix C is made up of the terms, Cij s.t.
% C44 = p*g*V*yb - m*g*yg
% C45 = -g*(p*V*xb - m*xg)
% C65 = -g*(p*V*zb - m*zg)
% C66 = p*g*V*yb = m*g*yg
% p is the water density
% g is the acceleration due to gravity
% V is the displaced volume
% center of gravity, cog is [xg,yg,zg]
% center of buoyancy, cob is [xb,yb,zb]
% m is the total glider's mass
%
% This is from a textbook where it is assumed that the origin is the center
%  of flotation which it very much is not as currently defined
    [N,~] = size(pm);    

    % Center of gravity
    CG = zeros(3,1); 
    for j = 1:3
        CG(j) = sum(pm(1:N,j+1)) / N;
    end

    % Total Volume
    V = volume(model.Mesh);
    [~,nn] = size(model.Mesh.Nodes);
    
    % Center of Buoyancy
    CB = sum(model.Mesh.Nodes,2)/nn;

    m = sum(pm(1:N,1));

    % This is only right if the center of buoyancy is at (0,0,0)
    % In other words it is not right
    C = zeros(6,6);
    C(2,2) = g*p*S;
    C(4,4) = g*(p*(S33 + V*CB(2)) - m*CG(2));
    C(4,5) = -g*(p*V*CB(1) - m*CG(1));
    C(6,5) = -g*(p*V*CB(3) - m*CG(3));
    C(6,6) = g*(p*(S11 + V*CB(2)) - m*CG(2));
end

