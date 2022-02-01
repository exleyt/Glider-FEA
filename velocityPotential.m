function [phi] = velocityPotential(S,C,N,N6,T,K)
%CALCULATE_VELOCITY_POTENTIAL_VECTOR Summary of this function goes here
%   Detailed explanation goes here

    % Defines the surface integral terms
    Gnks = zeros(S,S); % All N^2 Gnk
    Gnks_sum = zeros(S,6); % All N sums of Gnk * Fk as 6-dimensional vectors 
    Mnks = zeros(S,S); % All N^2 Mnk

    % Fills Gnks and Mnks
    parfor k = 1:S
        % Paremeterizes triangle to (u,v)
        Tk = T(:,:,k);
        ru = Tk(2,:) - Tk(1,:);
        rv = Tk(3,:) - Tk(1,:);

        % The area of the traingle
        A = 0.5*norm(cross(ru,rv));

        oneSixth = 1/6;
        
        % The guassian points to evaluate at
        r3 = [oneSixth*ru + oneSixth*rv; 
              2/3*ru + oneSixth*rv; 
              oneSixth*ru + 2/3*rv] + Tk(1,:);
        r7 = [0.1012865073235*ru + 0.1012865073235*rv; 
              0.7974269853531*ru + 0.1012865073235*rv;
              0.1012865073235*ru + 0.7974269853531*rv;
              0.4701420641051*ru + 0.0597158717898*rv;
              0.4701420641051*ru + 0.4701420641051*rv;
              0.0597158717898*ru + 0.4701420641051*rv;
              0.3333333333333*ru + 0.3333333333333*rv];
        
        % Guassian weights
        w3 = 1/3;
        w7 = [0.1259391805448;
            0.1259391805448;
            0.1259391805448;
            0.1323941527885;
            0.1323941527885;
            0.1323941527885;
            0.225];
        
        for n = 1:S 
            if n~=k
                % Sums each guassian function evaulation
                for m = 1:3
                    [fg,fm] = greenFunctionAndPartialXINormal(C(n,:),r3(m,:),N(k,:),K);
                    Gnks(n,k) = Gnks(n,k) + fg*w3;
                    Mnks(n,k) = Mnks(n,k) + fm*w3;
                end
            % More intensive guassian for the case when n=k
            else
                % Sums each guassian function evaulation
                for m = 1:7
                    [fg,fm] = greenFunctionAndPartialXINormal(C(n,:),r7(m,:),N(k,:),K);
                    Gnks(n,k) = Gnks(n,k) + fg*w7(m);
                    Mnks(n,k) = Mnks(n,k) + fm*w7(m);
                end
            end
            Gnks(n,k) = A*Gnks(n,k);
            Mnks(n,k) = A*Mnks(n,k);
        end
    end
    for n = 1:S
        for k = 1:S
            Gnks_sum(n,:) = Gnks_sum(n,:) + Gnks(n,k) * N6(k,:);
        end
        Mnks(n,n) = Mnks(n,n) + 2*pi;
    end
    
    phi = zeros(S,6);
    for j = 1:6
        phi(:,j) = linsolve(Mnks,Gnks_sum(:,j));
    end
end