function [fg,fm] = surfaceIntegralGMInner(xn,r,u,n,K,x,w)
%ESTIMATE_INNER_SURFACE_INTEGRAL_GM Summary of this function goes here
%   Detailed explanation goes here
    fg = 0;
    fm = 0;
    dvdx = (1 - u)*0.5;
    for i = 1:4
        [f,df] = greenFunctionAndPartialXINormal(xn, ...
            r(u,dvdx*x(i) + dvdx),n,K);
        fg = fg + f * w(i);
        fm = fm + df * w(i);
    end
    fg = fg*dvdx;
    fm = fm*dvdx;
end

