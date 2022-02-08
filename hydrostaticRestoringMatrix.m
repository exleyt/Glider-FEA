function [C] = hydrostaticRestoringMatrix(pm,g)
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

    cog = zeros(3,1); 
    for j = 1:3
        cog(j) = sum(pm(1:N,j+1)) / N;
    end

    m = sum(pm(1:N,1));

    % This is only right if the center of buoyancy is at (0,0,0)
    % In other words it is not right
    C = zeros(6,6);
    C(4,4) = -m*g*cog(2);
    C(4,5) = g*m*cog(1);
    C(6,5) = g*m*cog(3);
    C(6,6) = -m*g*cog(2);
end

