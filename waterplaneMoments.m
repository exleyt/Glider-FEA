function [S,CF,S11,S22] = waterplaneMoments(CL,P)
% Returns the waterplane area and moments of a list of triangles
%
% Given a list of point P(n,3) and connectivity list of triangles CL(N,3)
%  returns: 
% S the area of the waterplane
% CF the waterplane's center of flotation [S1/S,S2/S,0] where:
%  Sj = the surface integral of xj
% S11 and S22 where:
%  Sij =  the surface integral of xi*xj
    S = 0;
    CF = zeros(3,1); % [S1,S2,0]/S
    S11 = 0;
    S22 = 0;

    [N,~] = size(CL);
    for j = 1:N
        % Parameterized triangle
        r0 = P(CL(j,1),:);
        ru = P(CL(j,2),:) - r0;
        rv = P(CL(j,3),:) - r0;
        
        % Area of parallelogram
        A = norm(cross(ru,rv));

        S = S + 0.5*A;

        % S1
        CF(1) = CF(1) + A*(3*r0(1) + ru(1) + rv(1))/6;
        % S2
        CF(2) = CF(2) + A*(3*r0(2) + ru(2) + rv(2))/6;
        
        S11 = S11 + A*(6*r0(1)^2 + 4*(ru(1) + rv(1))*r0(1) + ru(1)^2 + ...
            ru(1)*rv(1) + rv(1)^2) / 12;
        S22 = S22 + A*(6*r0(2)^2 + 4*(ru(2) + rv(2))*r0(2) + ru(2)^2 + ...
            ru(2)*rv(2) + rv(2)^2) / 12;
    end

    CF = CF/S;
end

