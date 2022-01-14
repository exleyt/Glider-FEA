function [F] = calculate_exciting_forces_vector(T,phi,N,p,k,g,w,theta)
% Calculates the vector of exciting forces on the glider
%
% Force is the 6 dimensional vector defined as ??? 
% phi is the radiation potential vector (6 degrees of freesom)
% N is the normal vector
% p is the water density
% k is the wave number
% g is the acceleration due to gravity
% w is the gliders angular frequency
% theta is the direction of the incident wave
    F = zeros(6,1); % Exciting force vector
    [S,~] = size(phi);
    cost = cos(theta);
    sint = sin(theta);
    for m = 1:S
        ru = T(2,:,m) - T(1,:,m);
        rv = T(3,:,m) - T(1,:,m);
        n = norm(cross(ru,rv));
        A = k*(ru(3) - 1i*(ru(1)*cost - ru(2)*sint));
        B = k*(rv(3) - 1i*(rv(1)*cost - rv(2)*sint));
        C = A - B;
        expA = exp(A);
        expB = exp(B);
        expC = exp(C);
        D = k*(N(1)*cost + N(2)*sint + 1i*N(3));
        E = exp(-1i*k*(T(1,1,m)*cost + T(1,2,m)*sint + T(1,3,m))); % missing e^iwt term from my notes on purpose for now    
        G = (B*expA + A*expB - C)/(A*B*(-C));
        for j = 1:3
            F(j) = F(j) - 1i*p*g*n*E*(1i*N(j) - phi(j)*D)*G;
        end
        Hp = (B*expA + A*expB - C)/(A*C);
        Hu = (expC*(C - 1) - (expA*(A-1) + 1)/(A*A));
        Hv = ((B - 1)*(expC - 1)*expB/(B*C) ...
            - expB*(expC*(C - 1) + 1)/(C*C) + (expA - 1)/(A*B));
        for j = 4:6
            c = mod(j + 1,3) + 1;
            d = mod(j,3) + 1;
            Np = N(c)*T(1,c,m) - N(d)*T(1,d,m);
            Nu = N(c)*ru(d) - N(d)*ru(c);
            Nv = N(c)*rv(d) - N(d)*rv(c);
            F(j) = F(j) - 1i*w*p*n*E/B*((Np - phi(j)*D)*Hp + Nu*Hu + Nv*Hv);
        end
    end
end