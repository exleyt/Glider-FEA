function [Y] = surface_integral_of_green_function(xn,txi,X)
%SURFACE_INTEGRAL_OF_GREEN_FUNCTION Summary of this function goes here
%   Detailed explanation goes here
    addpath('asq')

    % Parametric function of triangular surface plane
    % Triangle defined on 0<=u<=1 and 0<=v<=1-u
    a = txi(2,:) - txi(1,:);
    b = txi(3,:) - txi(1,:);
    r = @(u,v) txi(1,:) + u*a + v*b;
    
    % Asq values (adjusted for testing speed)
    e = 10^-3;
    d = 20;
    
    % Surface integral integrand
    Gv = @(K,u,v) evaluate_green_function(xn,r(u,v),K,e,d);
    
    % Inner surface integral over v
    Gu = @(K,u) asq_3(Gv,K,u,0,1-u,e,d);
    
    % Outer surface integral over u
    Gnk = @(K) asq_2(Gu,K,0,1,e,d);
    
    % The nubmer of points to evaluate along K
    [~,s] = size(X); 
    
    % A matrix to contain each step of K and the surface integral at that value
    Y = zeros(s,1);
    
    % Fills that matrix with evaluations of Gnk at K
    for j = 1:s
        Y(j) = Gnk(X(j));
    end
end

