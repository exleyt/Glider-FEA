function [force] = calculate_exciting_forces_vector()
% Calculates the vector of exciting forces on the glider
%
% Force is the 6 dimensional vector defined as:
% F(j) = -iwp*Integral(n(j)*phi(j) - phi(i)*Partial(phi1 wrt n), S, dS)
% Where 
% w is the gliders angular velocity
% p is the water density
% n is the normal vector (6 degrees of freedom)
% phi is the radiation potential vector (6 degrees of freesom)
% phi1 is the velocity potential of incident wave with unit amplitude
    force = zeros(6,1);
end
% TODO
% Inputs
% sym ?
% partial derivative phi1 wrt n ?
% update summary
% write the function :P

