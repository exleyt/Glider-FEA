function [D] = addedMassAndDampingMatrices(Tri,phi,FN6,p)
% Calculates the added mass coefficient and damping coefficient matrices
% 
% The added mass coefficient matrix is the real component of V
% The damping coefficient matrix is an imaginary component of V * -w
% Vij = p*SurfaceIntegral(n(i)*phi(j),S,dS)
% Where:
% Tri is a (3,3,N) matrix of triangles where each row is [x,y,z]
% phi is a list of radiation potential vectors (6,N)
% FN is a list of normal vectors (S,3)
% p is the density of water at the operating temperature
% w is the gliders angular frequency
    [~,N] = size(phi);

    D = zeros(6,6);
    for k = 1:N
        % Parameterized triangle s.t. r = r0 + ru*u + rv*v
        r0 = Tri(1,:,k);
        ru = Tri(2,:,k) - r0;
        rv = Tri(3,:,k) - r0;

        % Area of traingle bounded by ru and rv
        A = 0.5*norm(cross(ru,rv));

        for j = 1:6
            Aphi = A*phi(j,k);

            % Simple integral for i:[1,3] s.t. n(i) = FN(k,i) 
            for i = 1:6
               D(i,j) = D(i,j) + Aphi*FN6(k,i);
            end
        end
    end

    D = p*D;
end