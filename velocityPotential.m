function [phi] = velocityPotential(CP,FN,Tri,K)
% Returns the time independent velocity potential across the surface
%
% Estimates the time independent velocity potential of a surface as a
%  piece-wise matrix (6,N) for each of the N triangles that make up the
%  surface at their center points using guassians to estimate surface 
%  integrals defined over green functions where:
% CP is a (N,3) matrix of triangle center points
% FN is a (N,3) matrix of triangle normals
% Tri is a (3,3,N) matrix of triangles where each row is [x,y,z]
% K = w^2/g
% w is the waves angular frequency
% g is acceleration due to gravity 
%
% Solves the linear equation:
% 2*pi*phi{j}(CP{n}) + Sum(k:[1,N],M{n,k}*phi{j}(CP{k})) = ...
%  Sum(k:[1,N],G{n,k}*FN6(CP{k})) s.t. n = 1...N
% G{nk} = surface integral of G(CP{n},X{k}) w.r.t Tri{k}
% M{nk} = surface integral of d(G(CP{n},X{k}))/d(FN{k}) w.r.t Tri{k}
    % N is the number of surface triangles
    [N,~] = size(CP);

    % Defines the surface integral terms
    Gnks = zeros(N,N); % All N^2 Gnk
    Gnks_sum = zeros(N,6); % All N sums of Gnk*FNk as 6-dimensional vectors 
    Mnks = zeros(N,N); % All N^2 Mnk

    % Fills Gnks and Mnks
    parfor k = 1:N
        % Paremeterizes triangle to (u,v) s.t. r = Tk(1,:) + ru*u + rv*v
        Tk = Tri(:,:,k);
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
        
        % Finds and stores the integral evaluation for each Gnk and Mnk
        for n = 1:N 
            if n~=k
                % Sums each guassian function evaulation
                for m = 1:3
                    [fg,fm] = greenFunctionAndPartialXINormal(CP(n,:),r3(m,:),FN(k,:),K);
                    Gnks(n,k) = Gnks(n,k) + fg*w3;
                    Mnks(n,k) = Mnks(n,k) + fm*w3;
                end
            % More intensive guassian for the case when n=k
            else
                % Sums each guassian function evaulation
                for m = 1:7
                    [fg,fm] = greenFunctionAndPartialXINormal(CP(n,:),r7(m,:),FN(k,:),K);
                    Gnks(n,k) = Gnks(n,k) + fg*w7(m);
                    Mnks(n,k) = Mnks(n,k) + fm*w7(m);
                end
            end
            Gnks(n,k) = A*Gnks(n,k);
            Mnks(n,k) = A*Mnks(n,k);
        end
    end

    % Creates each triangle's 6-dimensional normal vector
    FN6 = normal6DOF(CP,FN);
    
    for n = 1:N
        for k = 1:N
            % Sums Gnks*nk for each n at each dimension 
            Gnks_sum(n,:) = Gnks_sum(n,:) + Gnks(n,k) * FN6(k,:);
        end
        % Adds 2pi to each Mnk where k=n
        Mnks(n,n) = Mnks(n,n) + 2*pi;
    end
    
    % Solves the linear equation Sum(Mnk*phijk,k:[1,n]) = Gnks_sumnj for
    %  each dimesnion j
    phi = zeros(6,N);
    for j = 1:6
        phi(j,:) = linsolve(Mnks,Gnks_sum(:,j));
    end
end