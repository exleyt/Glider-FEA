addpath('asq');

% Creates a new pde object for mesh generation
model = createpde();

% Generates a discrete geometry object from stl file
importGeometry(model,"template-model\files\CORNER.STL");

% Moves the geometry underneath the waterline (Z=0)
model.Geometry.translate([0,0,-20]);

% Generates a finite element tetrahedra mesh object of the geometry
generateMesh(model, 'GeometricOrder','linear','Hmin',15);

% Plots mesh 
pdeplot3D(model);

% Makes a traingulation object of mesh
mto = triangulation(model.Mesh.Elements.', model.Mesh.Nodes.');

% Makes a traingle connectivity and point list of surface triangles
[L,P] = freeBoundary(mto);

% Makes a traingulation object of mesh surface
to = triangulation(L,P);

% Plots traingulation object
%figure();
%trimesh(to);
%hold on
%axis equal

% Sets S equal to the number of surface triangles
[S,~] = size(L);

% Finds each triangle's center and normal vector
C = incenter(to);
N = faceNormal(to);

% Plotes the triangulation object's normals
%quiver3(C(:,1),C(:,2),C(:,3), ...
%     F(:,1),F(:,2),F(:,3),0.5,'color','r');

% Finds each triangle's 6-dimensional normal vector
N6 = zeros(S,3);
N6(:,1:3) = N(:,:);
N6(:,4:6) = cross(C(:,:),N(:,:));

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

% Defines K from T=3:30 seconds
%K1 = (2*pi()/30)^2/9.8;
%K2 = (2*pi()/3)^2/9.8;

% The step of each estimation point of K
%x = 20;
%XS = (K2-K1)/(x-1); 
%X = K1:XS:K2;

% Defines a list of N triangles so that parfor can nicely distribute data
T = zeros(3,3,S);
for j = 1:S
    T(:,:,j) = to.Points(to.ConnectivityList(j,:),:);
end

phi = calculate_velocity_potential_vector(S,C,N,N6,T,K);

[A,B] = calculate_added_mass_and_damping_matrices(T,phi,N,p,w);
[F] = calculate_exciting_forces_vector(T,phi,N,p,k,g,w,theta); % missing exp(i*w*t) for now. Should be able to factor out term for use in solving eta
