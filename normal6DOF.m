function [FN6] = normal6DOF(CP,FN)
% Finds the 6 dof normal vectors of the discretized object
% 
% Finds the 6 dof normal vectors of the discretized object with respect to
%  that the center point of each discrete element where
% CP is the (N,3) list of center points
% FN is the (N,3) list of normal vectors
    [N,~] = size(CP);

    % Creates each triangle's 6-dimensional normal vector wrt the CP
    FN6 = zeros(N,6);
    FN6(:,1:3) = FN(:,:);
    FN6(:,4:6) = cross(CP(:,:),FN(:,:));
end

