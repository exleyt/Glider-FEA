function [a1,au,av] = surfIntPhiI(su,sv)
% Integral of [1,u,v]*e^(su*u)*e^(sv*v)dvdu s.t. u:[0,1] v:[0,1-u]
% 
% Since the anti-derivative of the integrand has terms 1/su, 1/sv, and 
%  1/(su-sv) the result must be split into five seperate cases to avoid 
%  1/0 errors.
% Alternative approaches include taking the limit of the anti-derivative
%  as the respective values approach su and sv or plugging su and sv into
%  the integrand before integration.
% However, this approach is by far the fastest. 
    if su == 0 && sv == 0
        % Limit of equ1 as (su,sv)->(0,0)
        a1 = 0.5;
        % Limit of equ2 as (su,sv)->(0,0)
        au = 1/6;
        % Limit of equ3 as (su,sv)->(0,0)
        av = 1/6;
    elseif su == 0
        % Limit of equ1 as su->0
        a1 = (exp(sv) - sv - 1)/(sv^2);
        % Limit of equ2 as su->0
        au = (2*exp(sv) - sv^2 - 2*sv - 2)/(2*sv^3);
        % Limit of equ3 as su->0
        av = ((sv - 2)*exp(sv) + sv + 2)/(sv^3);
    elseif sv == 0
        % Limit of equ1 as sv->0
        a1 = (exp(su) - su - 1)/(su^2);
        % Limit of equ2 as sv->0
        au = ((su - 2)*exp(su) + su + 2)/(su^3);
        % Limit of equ3 as sv->0
        av = (2*exp(su) - su^2 - 2*su - 2)/(2*su^3);
    elseif su == sv
        % Limit of equ1 as su->sv
        a1 = ((sv - 1)*exp(sv) + 1)/(sv^2);
        % Limit of equ2 as su->sv
        au = ((sv^2 - 2*sv + 2)*exp(sv) - 2)/(2*sv^3);
        % Limit of equ3 as su->sv
        av = ((sv^2 - 2*sv + 2)*exp(sv) - 2)/(2*sv^3);
    else
        % Integral of e^(su*u)*e^(sv*v)dvdu s.t. u:[0,1] v:[0,1-u]; equ1
        a1 = (exp(su)*sv - su*(exp(sv) - 1) - sv)/(su^2*(su - sv)); 
        % Integral of u*e^(su*u)*e^(sv*v)dvdu s.t. u:[0,1] v:[0,1-u]; equ2
        au = ((su^2 - su*(sv + 2) + sv)*exp(su)*sv + su^2*(exp(sv) - 1) + 2*su*sv - sv^2)/(su^2*(su-sv)^2*sv);
        % Integral of v*e^(su*u)*e^(sv*v)dvdu s.t. u:[0,1] v:[0,1-u]; equ3
        av = ((exp(su) - 1)*sv^2 - su^2*((sv - 1)*exp(sv) + 1) + su*sv*((sv - 2)*exp(sv) + 2))/(su*(su-sv)^2*sv^2);
    end
end