addpath('asq')
e = 10 ^ -6;
d = 20;
% Verticies of triangular surface
p1 = [0;0;0];
p2 = [1;0;-2];
p3 = [0;1;-2];
% Defines center point of velocity potential
p4 = [0.25;0.25;-2];
% Parametric function of triangular surface plane
% Triangle defined on 0<=u<=1 and 0<=v<=1-u
a = p2 - p1;
b = p3 - p1;
r = @(u,v) p1 + u*a + v*b;
% Surface integral integrand
Gv = @(K,u,v) evaluate_green_function(p4, r(u,v), K);
% Inner surface integral over v
Gu = @(K,u) asq_3(Gv,K,u,0,1-u,e,d);
% Outer surface integral over u
Mnk = @(K) asq_2(Gu,K,0,1,e,d);

