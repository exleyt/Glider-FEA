addpath('asq')

% Verticies of triangular surface
p1 = [19.9979904495083;23.6567196630356;-20];
p2 = [23.6567196630356;19.9979904495083;-20];
p3 = [22.0673361911350;22.0489968166775;-15.2187076286611];

n = [0.705464701512395,0.705464701512396,-0.0681110111513105];

% Defines center point of velocity potential
p4 = [21.9060437777891;21.8973319368885;-18.4601853230246];

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
Mv = @(K,u,v) evaluate_green_function_partial_xinormal(xn,r(u,v),K,e,d,n);

% Inner surface integral over v
Gu = @(K,u) asq_3(Gv,K,u,0,1-u,e,d);
Mu = @(K,u) asq_3(Mv,K,u,0,1-u,e,d);

% Outer surface integral over u
Gnk = @(K) asq_2(Gu,K,0,1,e,d);
Mnk = @(K) asq_2(Mu,K,0,1,e,d);

% Defines K from T=3:30 seconds
K1 = (2*pi()/30)^2/9.8;
K2 = (2*pi()/3)^2/9.8;

% One less than the nubmer of points to evaluate along K
s = 19;

% Displays approximate runtime
tG = timeit(@() Gnk(K2));
tM = timeit(@() Mnk(K2));
disp(['A single Gnk takes ~',sprintf('%.5f',tG),' seconds']);
disp(['A single Mnk takes ~',sprintf('%.5f',tM),' seconds']);
disp(['This process will take ~',sprintf('%.5f',(tG + tM)*(s+1)),' seconds']);

% The step of each point of K
KS = (K2-K1)/s;
KSL = K1:KS:K2;

% A matrix to contain each step of K and the surface integral at that value
GnK = zeros(s+1,1);
MnK = zeros(s+1,1);

% Fills that matrix with evaluations of Gnk at K
for j = 1:s+1
    GnK(j) = Gnk(KSL(j));
end

% Fills that matrix with evaluations of Mnk at K
for j = 1:s+1
    MnK(j) = Mnk(KSL(j));
end

% Plots all points to new window
figure()
hold on

% Plots the real and imaginary values of Gnk(K) with accuracy s
plot(KSL(:), real(GnK(:)), 'x');
plot(KSL(:), imag(GnK(:)), 'o');

% Plots all points to new window
figure()
hold on

% Plots the real and imaginary values of Mnk(K) with accuracy s
plot(KSL(:), real(MnK(:)), 'x');
plot(KSL(:), imag(MnK(:)), 'o');