function [Gnk,Mnk] = estimate_surface_integral_GM(xn,txi,n,K,x,w)
%ESTIMATE_SURFACE_INTEGRAL_GM Summary of this function goes here
%   Detailed explanation goes here

    % Parametric function of triangular surface plane
    % Triangle defined on 0<=u<=1 and 0<=v<=1-u
    ru = txi(2,:) - txi(1,:);
    rv = txi(3,:) - txi(1,:);
    r = @(u,v) txi(1,:) + u*ru + v*rv;
    m = norm(cross(ru,rv));
    
    Gnk = 0;
    Mnk = 0;
    for i = 1:4
        [fg,fm] = estimate_inner_surface_integral_GM(xn,r, ...
            x(i)*0.5 + 0.5,n,K,x,w);
        Gnk = Gnk + fg * w(i);
        Mnk = Mnk + fm * w(i);
    end
    Gnk = Gnk*0.5*m;
    Mnk = Mnk*0.5*m;
end

