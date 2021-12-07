function [Y] = surface_integral_of_green_function_partial_xinormal(xn,txi,X,n)
%SURFACE_INTEGRAL_OF_GREEN_FUNCTION_PARTIAL_XINORMAL_DERIVATIVE Summary of this function goes here
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
    Mv = @(K,u,v) evaluate_green_function_partial_xinormal(xn,r(u,v),K,e,d,n);
    
    % Inner surface integral over v
    Mu = @(K,u) asq_3(Mv,K,u,0,1-u,e,d);
    
    % Outer surface integral over u
    Mnk = @(K) asq_2(Mu,K,0,1,e,d);
    
    % The nubmer of points to evaluate along K
    [~,s] = size(X); 
    
    % A matrix to contain each step of K and the surface integral at that value
    Y = zeros(s,1);
    
    % Fills that matrix with evaluations of Mnk at K
    for j = 1:s
        Y(j) = Mnk(X(j));
    end
end

