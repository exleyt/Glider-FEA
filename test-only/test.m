addpath('asq');

model = createpde();
importGeometry(model,"template-model\files\CORNER.STL");
model.Geometry.translate([0,0,-20]);
generateMesh(model, 'GeometricOrder','linear','Hmin',20);
mto = triangulation(model.Mesh.Elements.', model.Mesh.Nodes.');
[T, P] = freeBoundary(mto);
to = triangulation(T, P);
[N,~] = size(T);
C = incenter(to);
F = faceNormal(to);
D = 6;
F6 = zeros(N,D);
F6(:,1:3) = F(:,:);
F6(:,4:6) = cross(C(:,:),F(:,:));
K1 = (2*pi()/30)^2/9.8;
K2 = (2*pi()/3)^2/9.8;
x = [sqrt(3/7 - 2/7 * sqrt(6/5)); -sqrt(3/7 - 2/7 * sqrt(6/5)); ...
     sqrt(3/7 + 2/7 * sqrt(6/5)); -sqrt(3/7 + 2/7 * sqrt(6/5))];
w = [(18 + sqrt(30))/36;(18 + sqrt(30))/36; ...
     (18 - sqrt(30))/36;(18 - sqrt(30))/36];
e = 10^-2;
e2 = 10^-3;
d = 20; 

k = 1;
n = 1;
txi = to.Points(to.ConnectivityList(k,:),:);
ru = txi(2,:) - txi(1,:);
rv = txi(3,:) - txi(1,:);
r = @(u,v) txi(1,:) + u*ru + v*rv;
m = norm(cross(ru,rv));
xn = C(n,:);
n = F(k,:);
u=0;

G = zeros(101);
M = zeros(101);

K = zeros(101);
