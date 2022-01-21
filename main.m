% NOTE: FFR -> Flagged For Removal
% NOTE: Replace this with iterative asq???
addpath('asq');

% Creates a new pde object for mesh generation
model = createpde();

% Generates a discrete geometry object from stl file
importGeometry(model,"template-model\files\CORNER.STL");

% Moves the geometry underneath the waterline (Z=0)
% NOTE: Make this done automatically???
model.Geometry.translate([0,0,-20]);

% Generates a finite element tetrahedra mesh object of the geometry
% NOTE: Hmin and Hmax should be set by user
generateMesh(model, 'GeometricOrder','linear','Hmin',15);

% FFR: Plots mesh
% pdeplot3D(model);

% Makes a traingulation object of the mesh
MTO = triangulation(model.Mesh.Elements.', model.Mesh.Nodes.');

% Makes a traingle connectivity and point list of surface triangles
[TC,P] = freeBoundary(MTO);

% Makes a traingulation object of surface triangles
to = triangulation(TC,P);

% FFR: Plots traingulation object
%figure();
%trimesh(to);
%hold on
%axis equal

% Records the number of surface triangles
[N,~] = size(to);

% Finds each triangle's center point and face normal vectors
CP  = incenter(to);
FN = faceNormal(to);

% FFR: Plotes the triangulation object's normals
%quiver3(C(:,1),C(:,2),C(:,3), ...
%     F(:,1),F(:,2),F(:,3),0.5,'color','r');

% Creates each triangle's 6-dimensional normal vector
FN6 = zeros(N,3);
FN6(:,1:3) = FN(:,:);
FN6(:,4:6) = cross(CP(:,:),FN(:,:));

% Plotes the triangulation object's higher dimension normals
%quiver3(C(:,1),C(:,2),C(:,3), ...
%     F6(:,4),F6(:,5),F6(:,6),0.5,'color','g');

% Setup new figure
%figure();
%hold on

% Temporarily define response inputs
w = 2*pi()/30;
theta = 0;
g = 9.8;
p = 1; 
K = w^2/g;
k = 1; 
a = 1;
pm = [5,5;6,0;8,0;10,0];

% FFR: Defines K from T=3:30 seconds
%K1 = (2*pi()/30)^2/9.8;
%K2 = (2*pi()/3)^2/9.8;

% FFR: The step of each estimation point of K
%x = 20;
%XS = (K2-K1)/(x-1); 
%X = K1:XS:K2;

% Defines a list of N triangles so that parfor can nicely distribute data
T = zeros(3,3,N);
for j = 1:N
    T(:,:,j) = to.Points(to.ConnectivityList(j,:),:);
end

phi = velocityPotential(N,CP,FN,FN6,T,K);
[A,B] = AddedMassAndDampingMatrices(T,phi,FN,p,w);
[F] = excitingForce(T,phi,FN,p,k,g,w,theta);
M = bodyInertiaMatrix(pm);
C = hydrostaticRestoringMatrix(pm,g);

nu = linsolve(-w^2*(M + A) + 1i*w*B + C,F);
H = nu / a;