function [f,df] = greenFunctionAndPartialXINormal(x,xi,FN,K)
% Returns the green function and its partial normal of xi
%
% Finds the complex valued influence of a pulsating source located at xi on
%  the potential at x and its gradient of xi dotted with the normal where:
% x is the location of the potential
% xi is the location of the pulsating source
% FN is the normal vector of the tile corresponding to xi 
% K = w^2/g
% w is the waves angular frequency
% g is acceleration due to gravity
    f = greenFunction(x,xi,K); % the green function evaluated at xi
    
    % Finds the xi gradient of f
    grad = zeros(3,1);
    for j = 1:3
        grad(j) = greenFunctionPartialXI(x,xi,K,j,f);
    end

    df = dot(real(grad), FN) + 1i*dot(imag(grad), FN);
end

