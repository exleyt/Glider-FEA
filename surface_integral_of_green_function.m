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
    Mnk = @(K) asq_2(Gu,K,0,1,e,d);
    
    % The nubmer of points to evaluate along K
    [~,s] = size(X); 
    s4 = s / 4; % 4* less resolution than our later makima estimation

    % The step of each point of K
    KS = (X(end)-X(1))/(s4-1);
    KSL = X(1):KS:X(end);
    
    % A matrix to contain each step of K and the surface integral at that value
    KMnk = zeros(s4,2);
    
    % Fills the matrix with values of K
    KMnk(1:s4,1) = KSL;
    
    % Fills that matrix with evaluations of Mnk at K
    for j = 1:s4
        KMnk(j,2) = Mnk(KMnk(j,1));
    end

    % Records the makima estimations of Mnk(K)
    Y = zeros(2,s);
    Y(1,:) = makima(KMnk(:,1), real(KMnk(:,2)), X);
    Y(2,:) = makima(KMnk(:,1), imag(KMnk(:,2)), X);
end

