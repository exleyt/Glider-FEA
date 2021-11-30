addpath('asq')
e = 10 ^ -6;
d = 20;
p1 = [0;0;0];
p2 = [1;0;-2];
p3 = [0;1;-2];
p4 = [0.25;0.25;-2];
a = p2 - p1;
b = p3 - p1;
r = @(u,v) p1 + u*a + v*b;
Gv = @(K,u,v) evaluate_green_function(p4, r(u,v), K);
Gu = @(K,u) asq_3(Gv,K,u,0,1-u,e,d);
Mnk = @(K) asq_2(Gu,K,0,1,e,d);
Mnk(0.04)
Mnk(0.22)
Mnk(0.4)