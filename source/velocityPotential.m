function [phi] = velocityPotential(Tri,CP,FN,K)
% Returns the time independent velocity potential across the surface
%
% Estimates the time independent velocity potential of a surface as a
%  piece-wise matrix (6,N) for each of the N triangles that make up the
%  surface at their center points using guassians to estimate surface 
%  integrals defined over green functions where:
% Tri is a (3,3,N) matrix of triangles where each row is [x,y,z]
% CP is a (N,3) matrix of triangle center points
% FN is a (N,3) matrix of triangle normals
% K = w^2/g
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
        r = [oneSixth*ru + oneSixth*rv; 
             4*oneSixth*ru + oneSixth*rv; 
             oneSixth*ru + 4*oneSixth*rv] + Tk(1,:);
        
        % Guassian weights
        w = 2*oneSixth;
        
        % Finds and stores the integral evaluation for each Gnk and Mnk
        for n = 1:N 
            if n ~= k
                % Sums each guassian function evaulation
                for m = 1:3
                    [fg,fm] = greenFunctionAndPartialXINormal(CP(n,:),r(m,:),FN(k,:),K);
                    Gnks(n,k) = Gnks(n,k) + fg*w;
                    Mnks(n,k) = Mnks(n,k) + fm*w;
                end
                Gnks(n,k) = A*Gnks(n,k);
                Mnks(n,k) = A*Mnks(n,k);
            else
                Mnks(n,k) = 2*pi;
            end
        end
    end

    % Creates each triangle's 6-dimensional normal vector
    FN6 = normal6DOF(CP,FN);
    
    for n = 1:N
        for k = 1:N
            % Sums Gnks*nk for each n at each dimension 
            Gnks_sum(n,:) = Gnks_sum(n,:) + Gnks(n,k) * FN6(k,:);
        end
    end
    
    % Solves the linear equation Sum(Mnk*phijk,k:[1,n]) = Gnks_sumnj for
    %  each dimesnion j
    phi = zeros(6,N);
    for j = 1:6
        phi(j,:) = linsolve(Mnks,Gnks_sum(:,j));
    end
end