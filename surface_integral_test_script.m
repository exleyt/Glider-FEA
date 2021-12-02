addpath('asq')

% Verticies of triangular surface
p1 = [16.9830640688435;9.84160996159739;-20];
p2 = [12.4200889907704;9.23805040555441;-20];
p3 = [15.1945555560609;14.3413934785717;-20];

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

% Inner surface integral over v
Gu = @(K,u) asq_3(Gv,K,u,0,1-u,e,d);

% Outer surface integral over u
Mnk = @(K) asq_2(Gu,K,0,1,e,d);

% Defines K from T=3:30 seconds
K1 = (2*pi()/30)^2/9.8;
K2 = (2*pi()/3)^2/9.8;

% The nubmer of points to evaluate along K
s = 40;
s4 = 4*s;

% Displays approximate runtime
t = timeit(@() Mnk(K2));
disp(['A single surface integral takes ~',sprintf('%.5f',t),' seconds']);
disp(['This process will take ~',sprintf('%.5f',t*s/60),' minutes']);

% The step of each point of K
KS = (K2-K1)/s;
KSL = K1:KS:K2;
KS4 = (K2-K1)/s4;
KS4L = K1:KS4:K2; % Simulate 4 times more accurate s

% A matrix to contain each step of K and the surface integral at that value
KMnk = zeros(s+1,1);

% Fills that matrix with evaluations of Mnk at K
for j = 1:s+1
    KMnk(j,1) = Mnk(KSL(j));
end

% Plots all points to new window
figure()
hold on

% Plots the real and imaginary values of Mnk(K) with accuracy s
plot(KSL(:), real(KMnk(:,1)), 'o');
plot(KSL(:), imag(KMnk(:,1)), 'o');

% Plots the real values of makima estimate of Mnk(K) with accuracy 4s
realy = makima(KSL(:), real(KMnk(:,1)), KS4L);
plot(KS4L(:),realy,"*");

% Plots the imaginary values of makima estimate of Mnk(K) with accuracy 4s
imagy = makima(KSL(:), imag(KMnk(:,1)), KS4L);
plot(KS4L(:),imagy,"*");

% A matrix to contain each step of K and the surface integral at that value
% This time with accuracy s4
KMnk = zeros(s4+1, 1);

% Fills that matrix with evaluations of Mnk at K
for j = 1:s4+1
    KMnk(j,1) = Mnk(KS4L(j));
end

% Plots the real and imaginary values of Mnk(K) with accuracy s4
plot(KS4L(:), real(KMnk(:,1)), '-');
plot(KS4L(:), imag(KMnk(:,1)), '-');