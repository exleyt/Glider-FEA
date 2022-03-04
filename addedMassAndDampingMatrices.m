function [A,B] = addedMassAndDampingMatrices(Tri,phi,FN,p,w)
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

    V = zeros(6,6);
    for k = 1:N
        % Parameterized triangle s.t. r = r0 + ru*u + rv*v
        r0 = Tri(1,:,k);
        ru = Tri(2,:,k) - r0;
        rv = Tri(3,:,k) - r0;

        % Area of traingle bounded by ru and rv * p
        Ap = 0.5*norm(cross(ru,rv))*p;

        for j = 1:6
            Apphi = Ap*phi(j,k);

            % Simple integral for i:[1,3] s.t. n(i) = FN(k,i) 
            for i = 1:3
               V(i,j) = V(i,j) + Apphi*FN(k,i);
            end

            % Less simple integral for i:[4,6] s.t. n(i) = n0 + nu*u + nv*v 
            for i = 4:6
                a = mod(i + 1,3) + 1;   % 3,1,2
                b = mod(i,3) + 1;       % 2,3,1

                 % n0 + nu*u + nv*v
                n0 = FN(k,a)*r0(b) - FN(k,b)*r0(a);
                nu = FN(k,a)*ru(b) - FN(k,b)*ru(a);
                nv = FN(k,a)*rv(b) - FN(k,b)*rv(a);
                                        
                V(i,j) = V(i,j) + Apphi*(n0 + (nu + nv)/3);
            end
        end
    end
    A = real(V);
    B = -imag(V)*w;
end