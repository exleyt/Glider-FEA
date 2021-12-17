addpath('asq');

% Creates a new pde object for mesh generation
model = createpde();

% Generates a discrete geometry object from stl file
importGeometry(model,"template-model\files\CORNER.STL");

% Moves the geometry underneath the waterline (Z=0)
model.Geometry.translate([0,0,-20]);

% Generates a finite element tetrahedra mesh object of the geometry
generateMesh(model, 'GeometricOrder','linear','Hmin',12);

% Plots mesh 
%pdeplot3D(model);

% Makes a traingulation object of mesh
mto = triangulation(model.Mesh.Elements.', model.Mesh.Nodes.');

% Makes a traingle connectivity and point list of surface triangles
[T, P] = freeBoundary(mto);

% Makes a traingulation object of mesh surface
to = triangulation(T, P);

% Plots traingulation object
%figure();
%trimesh(to);
%hold on
%axis equal

% Sets N equal to the number of surface triangles
[N,~] = size(T);

% Finds each triangle's center and normal vector
C = incenter(to);
F = faceNormal(to);

% Plotes the triangulation object's normals
%quiver3(C(:,1),C(:,2),C(:,3), ...
%     F(:,1),F(:,2),F(:,3),0.5,'color','r');

% Finds each triangle's 6-dimensional normal vector
D = 6; % Degrees of freedom
F6 = zeros(N,D);
F6(:,1:3) = F(:,:);
F6(:,4:6) = cross(C(:,:),F(:,:));

% Plotes the triangulation object's higher dimension normals
%quiver3(C(:,1),C(:,2),C(:,3), ...
%     F6(:,4),F6(:,5),F6(:,6),0.5,'color','g');

% Setup new figure
%figure();
%hold on

% Defines K from T=3:30 seconds
K1 = (2*pi()/30)^2/9.8;
K2 = (2*pi()/3)^2/9.8;

% The step of each estimation point of K
%x = 20;
%XS = (K2-K1)/(x-1); 
%X = K1:XS:K2;

% Defines the surface integral terms
Gnks = zeros(N,N); % All N^2 Gnk
Gnks_sum = zeros(N,D); % All N sums of Gnk * Fk as 6-dimensional vectors 
Mnks = zeros(N,N); % All N^2 Mnk

% Defines values for guasian quadrature
x = [sqrt(3/7 - 2/7 * sqrt(6/5)); -sqrt(3/7 - 2/7 * sqrt(6/5)); ...
     sqrt(3/7 + 2/7 * sqrt(6/5)); -sqrt(3/7 + 2/7 * sqrt(6/5))];
w = [(18 + sqrt(30))/36;(18 + sqrt(300))/36; ...
     (18 - sqrt(30))/36;(18 - sqrt(3))/36];

% Defines a list of N triangles
txk = zeros(3,3,N);
for k = 1:N
    txk(:,:,k) = to.Points(to.ConnectivityList(k,:),:);
end

tic
% Fills Gnks and Mnks
parfor k = 1:N
    for n = 1:N
        [Gnks(n,k),Mnks(n,k)] = estimate_surface_integral_GM(C(n,:), ...
            txk(:,:,k),F(k,:),K1,x,w);
    end
end
for n = 1:N
    for k = 1:N
        Gnks_sum(n,:) = Gnks_sum(n,:) + Gnks(n,k) * F6(k,:);
        Mnks(n,1) = Mnks(n,1) + 2*pi;
    end
end
toc
