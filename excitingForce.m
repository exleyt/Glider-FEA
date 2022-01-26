function [F] = excitingForce(T,phi,N,p,k,g,theta)
% Calculates the vector of exciting forces on the glider
%
% T is a list of triangles (3,3,S) 
% phi is a list of radiation potential vectors (S,6) 
% N is a list of normal vectors (S,3)
% p is the water density
% k is the wave number
% g is the acceleration due to gravity
% theta is the direction of the incident wave
    F = zeros(6,1);
    [S,~] = size(phi);
    cost = cos(theta);
    sint = sin(theta);
    for m = 1:S
        r0 = T(1,:,m);
        ru = T(2,:,m) - r0;
        rv = T(3,:,m) - r0;
        n = norm(cross(ru,rv));
        s1 = k*(r0(3) - 1i*(r0(1)*cost + r0(2)*sint));
        su = k*(ru(3) - 1i*(ru(1)*cost + ru(2)*sint));
        sv = k*(rv(3) - 1i*(rv(1)*cost + rv(2)*sint));
        exps1 = exp(s1);
        [a1,au,av] = surfIntPhiI(su,sv);
        sF1 = n*a1*exps1;
        SphiI = k*(N(m,1)*cost + N(m,2)*sint + 1i*N(m,3));
        for j = 1:3
            F(j) = F(j) + sF1*(1i*N(m,j) - phi(m,j)*SphiI);
        end
        sF4 = 1i*n*exps1;
        for j = 4:6
            a = mod(j + 1,3) + 1; % 3,1,2
            b = mod(j,3) + 1; % 2,3,1
            n0 = N(m,a)*r0(b) - N(m,b)*r0(a);
            nu = N(m,a)*ru(b) - N(m,b)*ru(a);
            nv = N(m,a)*rv(b) - N(m,b)*rv(a);
            F(j) = F(j) + sF4*(n0*a1 + nu*au + nv*av) - sF1*phi(m,j)*SphiI;
        end
    end
    F = -F*1i*p*g; % This is *e^iwt but I think it disappears when 
                   %  converting to frequency domain.
end