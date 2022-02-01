function [F] = excitingForce(T,phi,N,p,k,g,theta)
% Calculates the vector of exciting forces on the glider
%
% Calculates the surface integral of -i*w*p*(n(j)*phiI - phi(j)*dphiI/dn)
% T is a list of triangles (3,3,S) s.t. T(2,:,1) is the postiion vector 
%  [x,y,z] of the first triangle's second point 
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
        % Parameterized triangle s.t. r = r0 + ru*u + rv*v
        r0 = T(1,:,m);
        ru = T(2,:,m) - r0;
        rv = T(3,:,m) - r0;
        
        % Area of parallelogram bounded by ru and rv 
        n = norm(cross(ru,rv));

        % e^s1
        s1 = k*(r0(3) - 1i*(r0(1)*cost + r0(2)*sint));
        % e^su*u
        su = k*(ru(3) - 1i*(ru(1)*cost + ru(2)*sint));
        % e^sv*v
        sv = k*(rv(3) - 1i*(rv(1)*cost + rv(2)*sint));

        exps1 = exp(s1);
        
        % Integral of [1,u,v]*e^(su*u)*e^(sv*v)dvdu s.t. u:[0,1] v:[0,1-u]
        [a1,au,av] = surfIntPhiI(su,sv);  

        % Parts that factor out of the integral
        sF = n*exps1;
        % A part of dphiI/dn that is independed of (u,v) and doesn't factor
        SphiI = k*(N(m,1)*cost + N(m,2)*sint + 1i*N(m,3));

        % Facotrs out of whole integral for j:[1,3]
        sF13 = sF*a1;
        
        % Simple integral where for j:[1,3] s.t. dphiI/dn = N(m,j)
        for j = 1:3
            F(j) = F(j) + sF13*(1i*N(m,j) - phi(m,j)*SphiI);
        end

        % Facotrs out of whole integral for j:[4,6]
        sF46 = 1i*sF;

        % Less simple integral for j:[4,6] s.t. dphiI/dn = n0 + nu*u + nv*v
        for j = 4:6
            a = mod(j + 1,3) + 1; % 3,1,2
            b = mod(j,3) + 1; % 2,3,1
            n0 = N(m,a)*r0(b) - N(m,b)*r0(a);
            nu = N(m,a)*ru(b) - N(m,b)*ru(a);
            nv = N(m,a)*rv(b) - N(m,b)*rv(a);
            F(j) = F(j) + sF46*(n0*a1 + nu*au + nv*av) - sF13*phi(m,j)*SphiI;
        end
    end
    
    % This is F/e^iwt but I think e^iwt disappears in frequency domain.
    F = -F*1i*p*g; 
end