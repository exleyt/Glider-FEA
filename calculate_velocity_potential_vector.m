function [phi] = calculate_velocity_potential_vector(N,C,F,F6,T,K)
%CALCULATE_VELOCITY_POTENTIAL_VECTOR Summary of this function goes here
%   Detailed explanation goes here

    % Defines the surface integral terms
    Gnks = zeros(N,N); % All N^2 Gnk
    Gnks_sum = zeros(N,6); % All N sums of Gnk * Fk as 6-dimensional vectors 
    Mnks = zeros(N,N); % All N^2 Mnk
    
    % Defines values for guasian quadrature
    x = [sqrt(3/7 - 2/7 * sqrt(6/5)); -sqrt(3/7 - 2/7 * sqrt(6/5)); ...
         sqrt(3/7 + 2/7 * sqrt(6/5)); -sqrt(3/7 + 2/7 * sqrt(6/5))];
    w = [(18 + sqrt(30))/36;(18 + sqrt(300))/36; ...
         (18 - sqrt(30))/36;(18 - sqrt(3))/36];

    % Fills Gnks and Mnks
    parfor k = 1:N
        for n = 1:N
            [Gnks(n,k),Mnks(n,k)] = estimate_surface_integral_GM(C(n,:), ...
                T(:,:,k),F(k,:),K,x,w);
        end
    end
    for n = 1:N
        for k = 1:N
            Gnks_sum(n,:) = Gnks_sum(n,:) + Gnks(n,k) * F6(k,:);
        end
        Mnks(n,n) = Mnks(n,n) + 2*pi;
    end
    
    phi = zeros(N,6);
    for j = 1:6
        phi(:,j) = linsolve(Mnks,Gnks_sum(:,j));
    end
end

