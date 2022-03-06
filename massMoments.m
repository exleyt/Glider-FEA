function [m,Cg] = massMoments(pm)
% Returns the mass and center of mass of a list of point masses
%
% Where
% pm is the matrix of [mjs,pjs]
% mj is the jth point mass
% pj is the position vector of the jth point mass (xj,yj,zj)
% Returns
% m is the total mass
% Cg is the center of gravity
    m = sum(pm(:,1)); % 
    Cg = sum(pm(:,2:4).*pm(:,1),1)/m; % 
end

