function [amcMatrix,dcMatrix] = calculate_added_mass_and_damping_coefficient_matrices(w,p,n,phi)
% Summary
% 
% The added mass coefficient matrix is the real component of V
% The damping coefficient matrix is the imaginary component of V
% V = p*Integral(n(i)*phi(j), S, dS)
% Where:
% w is the gliders angular velocity
% p is the water density
% n is the normal vector (6 degrees of freedom)
% phi is the radiation potential vector (6 degrees of freesom)
    valueMatrix = zeros(6);
    for i = 1:6
        for j = 1:6
            valueMatrix(i,j) = p*n(i)*phi(j); % NOT COMPLETE
        end
    end  
    amcMatrix = real(valueMatrix);
    dcMatrix = -imag(valueMatrix)/w;
end
% TODO
% Take actual surface integral (over finite element mesh)
% Add mesh parameter
% sym w ?
