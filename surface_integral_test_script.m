addpath('asq')

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

% Asq values
e = 10^-3;
d = 20;

% Surface integral integrand
Gv = @(K,u,v) evaluate_green_function(p4,r(u,v),K,e,d);

% Inner surface integral over v
Gu = @(K,u) asq_3(Gv,K,u,0,1-u,e,d);

% Outer surface integral over u
Mnk = @(K) asq_2(Gu,K,0,1,e,d);

% Defines K from T=3:30 seconds
K1 = (2*pi()/30)^2/9.8;
K2 = (2*pi()/3)^2/9.8;

% The nubmer of points to evaluate along K
s = 20;

% Displays approximate runtime
t = timeit(@() Mnk(K2));
disp(['A single surface integral takes ~',sprintf('%.5f',t),' seconds']);
disp(['This process will take ~',sprintf('%.5f',t*s/60),' minutes']);

% The step of each point of K
KS = (K2-K1)/(s-1);

% A matrix to contain each step of K and the surface integral at that value
KMnk = zeros(s,2);

% Fills the matrix with values of K
KMnk(1:s,1) = K1:KS:K2;

% Fills that matrix with evaluations of Mnk at K
for j = 1:s
    KMnk(j,2) = Mnk(KMnk(j,1));
end

% Plots the imaginary values of Mnk(K)
hold on
for j = 1:s
    plot(KMnk(j,1),imag(KMnk(j,2)),"o")
end

% Plots the real values of Mnk(K)
figure()
hold on
for j = 1:s
    plot(KMnk(j,1),real(KMnk(j,2)),"x")
end
