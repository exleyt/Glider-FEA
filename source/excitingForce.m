function [F] = excitingForce(T,phi,FN6,p,K,g,theta)
% Calculates the vector of exciting forces on the glider
%
% Calculates the surface integral of -i*w*p*(n(j)*phiI - phi(j)*dphiI/dn)
% T is a list of triangles (3,3,S) s.t. T(2,:,1) is the postiion vector 
%  [x,y,z] of the first triangle's second point 
% phi is a list of radiation potential vectors (6,S) 
% FN is a list of normal vectors (S,3)
% p is the water density
% k is the wave number
% g is the acceleration due to gravity
% theta is the direction of the incident wave
    F = zeros(6,1);
    [~,S] = size(phi);
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
        s0 = K*(r0(3) - 1i*(r0(1)*cost + r0(2)*sint));
        % e^su*u
        su = K*(ru(3) - 1i*(ru(1)*cost + ru(2)*sint));
        % e^sv*v
        sv = K*(rv(3) - 1i*(rv(1)*cost + rv(2)*sint));

        exps1 = exp(s0);
        
        % Integral of [1,u,v]*e^(su*u)*e^(sv*v)dvdu s.t. u:[0,1] v:[0,1-u]
        [a1,~,~] = surfIntPhiI(su,sv);  

        % Parts that factor out of the integral
        sF = n*exps1*a1;

        % A part of dphiI/dn that is independed of (u,v) and doesn't factor
        SphiI = K*(FN6(m,1)*cost + FN6(m,2)*sint + 1i*FN6(m,3));
        
        % Simple integral for j:[1,3] s.t. n(j) = N(m,j)
        for j = 1:6
            F(j) = F(j) + sF*(1i*FN6(m,j) - phi(j,m)*SphiI);
        end
    end
    
    % This is F/e^iwt but I think e^iwt disappears in frequency domain.
    F = -F*1i*p*g; 
end