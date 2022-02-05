function [C] = hydrostaticRestoringMatrix(pm,g)
% A 6x6 matrix describing the body's hydrostatic restoring force
    C = zeros(6,6);
    [N,~] = size(pm);     
    cog = zeros(3,1); 
    for j = 1:3
        cog(j) = sum(pm(1:N,j+1)) / N;
    end
    m = sum(pm(1:N,1));
    C(4,4) = -m*g*cog(2);
    C(4,5) = g*m*cog(1);
    C(6,5) = g*m*cog(3);
    C(6,6) = -m*g*cog(2);
end

