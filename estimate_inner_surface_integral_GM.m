function [fg,fm] = estimate_inner_surface_integral_GM(xn,r,u,n,K,x,w,e,d)
%ESTIMATE_INNER_SURFACE_INTEGRAL_GM Summary of this function goes here
%   Detailed explanation goes here
    fg = 0;
    fm = 0;
    dvdx = (1 - u)*0.5;
    for i = 1:4
        [f,df] = evaluate_green_function_and_partial_xinormal(xn, ...
            r(u,dvdx*x(i) + dvdx),n,K,e,d);
        fg = fg + f * w(i);
        fm = fm + df * w(i);
    end
    fg = fg*dvdx;
    fm = fm*dvdx;
end

