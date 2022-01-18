function [A,B] = calculate_added_mass_and_damping_matrices(T,phi,N,p,w)
% Calculates the added mass coefficient and damping coefficient matrices
% 
% The added mass coefficient matrix is the real component of V
% The damping coefficient matrix is an imaginary scalar component of V
% V = p*SurfaceIntegral(N6(i)*phi(j),S,dS)
% Where:
% phi is the radiation potential vector
% N6 is the normal vector
% p is the density of water at the operating temperature
% w is the gliders angular frequency
    rhs = zeros(6,6);
    [S,~] = size(phi);
    for k = 1:S
        ru = T(2,:,k) - T(1,:,k);
        rv = T(3,:,k) - T(1,:,k);
        mp = p * norm(cross(ru,rv));
        r = T(1,:,k)*0.5 + ru/6 -rv/6;
        for j = 1:6
            for i = 1:3
               rhs(i,:) = rhs(i,:) + mp*phi(k,i)*r(i);
            end
            for i = 4:6
                c = mod(i + 1,3) + 1;
                d = mod(i,3) + 1;
                rhs(i,:) = rhs(i,:) + mp*phi(k,i)*(N(c)*r(d) ...
                    - N(d)*r(c));
            end
        end
    end
    A = real(rhs);
    B = -imag(rhs)*w;
end