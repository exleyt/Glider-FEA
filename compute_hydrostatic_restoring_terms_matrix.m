function [C] = compute_hydrostatic_restoring_terms_matrix(pm,g)
% A 6x6 matrix describing the body's hydrostatic restoring force
%
    C = zeros(6,6);
    cog = zeros(3,1); 
    for j = 1:3
        cog(j) = sum(pm(j+1,1:N)) / N;
    end
    m = sum(pm(1,1:N));
    C(4,4) = -m*g*cog(2);
    C(4,5) = g*m*cog(1);
    C(6,5) = g*m*cog(3);
    C(6,6) = -m*g*cog(2);
end

