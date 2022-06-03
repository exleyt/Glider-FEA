function [S,Sx,Sy,Sxy,Sxx,Syy] = waterplaneMoments(CL,P,FN)
% Returns the waterplane area and moments of a list of triangles
%
% Given a list of point P(n,3), connectivity list of triangles CL(N,3), and
%  list of normals FN(n,3) returns: 
% S the area of the waterplane
% Sx the integral of x over the waterplane area
% Sy the integral of y over the waterplane area
% Sxy the integral of xy over the waterplane area
% Sxx the integral of xx over the waterplane area
% Syy the integral of yy over the waterplane area
    S = 0;
    Sx = 0;
    Sy = 0;
    Sxy = 0;
    Sxx = 0;
    Syy = 0;
    [N,~] = size(CL);

    % Sums discretize elements of the waterplane
    for j = 1:N
        % Parameterized triangle
        r0 = P(CL(j,1),:);
        ru = P(CL(j,2),:) - r0;
        rv = P(CL(j,3),:) - r0;
        
        % Area of parallelogram
        A = norm(cross(ru,rv));

        S = S + 0.5*A*FN(j,3);
        Sx = Sx + A*FN(j,3)*(3*r0(1) + ru(1) + rv(1))/6;
        Sy = Sy + A*FN(j,3)*(3*r0(2) + ru(2) + rv(2))/6;
        Sxy = Sxy + A*FN(j,3)*(4*r0(1)*(3*r0(2) + ru(2) + rv(2)) ...
            + 4*r0(2)*(ru(1) + rv(1)) + ru(1)*(2*ru(2) + rv(2)) ...
            + rv(1)*(ru(2) + 2*rv(2)))/24;
        Sxx = Sxx + A*FN(j,3)*(6*r0(1)^2 + 4*(ru(1) + rv(1))*r0(1) + ru(1)^2 + ...
            ru(1)*rv(1) + rv(1)^2) / 12;
        Syy = Syy + A*FN(j,3)*(6*r0(2)^2 + 4*(ru(2) + rv(2))*r0(2) + ru(2)^2 + ...
            ru(2)*rv(2) + rv(2)^2) / 12;
    end
end

