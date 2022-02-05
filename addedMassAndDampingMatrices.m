function [A,B] = addedMassAndDampingMatrices(Tri,phi,FN,p,w)
% Calculates the added mass coefficient and damping coefficient matrices
% 
% The added mass coefficient matrix is the real component of V
% The damping coefficient matrix is an imaginary scalar component of V
% V = p*SurfaceIntegral(N6(i)*phi(j),S,dS)
% Where:
% phi is the radiation potential vector
% N is the normal vector
% p is the density of water at the operating temperature
% w is the gliders angular frequency
    rhs = zeros(6,6);
    [~,S] = size(phi);
    for k = 1:S
        r0 = Tri(1,:,k);
        ru = Tri(2,:,k) - r0;
        rv = Tri(3,:,k) - r0;
        Ap = 0.5*norm(cross(ru,rv))*p;
        for j = 1:6
            Apphi = Ap*phi(j,k);
            for i = 1:3
               rhs(i,j) = rhs(i,j) + Apphi*FN(k,i);
            end
            for i = 4:6
                a = mod(j + 1,3) + 1; % 3,1,2
                b = mod(j,3) + 1; % 2,3,1
                n0 = FN(k,a)*r0(b) - FN(k,b)*r0(a);
                nu = FN(k,a)*ru(b) - FN(k,b)*ru(a);
                nv = FN(k,a)*rv(b) - FN(k,b)*rv(a);
                rhs(i,j) = rhs(i,j) + Apphi*(n0 + (nu + nv)/3);
            end
        end
    end
    A = real(rhs);
    B = -imag(rhs)*w;
end