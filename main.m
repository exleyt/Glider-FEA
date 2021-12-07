% Creates a new pde object
model = createpde();

% Generates a discrete geometry object from stl file
importGeometry(model,"template-model\files\CORNER.STL");

% Moves the geometry underneath the waterline (Z=0)
model.Geometry.translate([0,0,-20]);

% Generates a finite element tetrahedra mesh object of the geometry
generateMesh(model, 'GeometricOrder','linear','Hmin',5);

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

% Declares the pde as a system of equations N equations
model.PDESystemSize = N;

% Finds each triangle's center and normal vector
C = incenter(to);
F = faceNormal(to);

% Plotes the triangulation object's normals
%quiver3(C(:,1),C(:,2),C(:,3), ...
%     F(:,1),F(:,2),F(:,3),0.5,'color','r');

% Setup new figure
%figure();
%hold on

% Defines K from T=3:30 seconds
K1 = (2*pi()/30)^2/9.8;
K2 = (2*pi()/3)^2/9.8;

% The step of each estimation point of K
x = 20;
XS = (K2-K1)/(x-1); 
X = K1:XS:K2;

for n = 1:1 % Systems of equations building loop
    Gnks = zeros(x,N);
    Mnks = zeros(x,N,N);
    for k = 1:N
        [Gnk] = surface_integral_of_green_function(C(n,:), to.Points(to.ConnectivityList(k,:),:),X);
        Gnks(:,n) = Gnks(:,n) + Gnk(:);
        [Mnks(:,n,k)] = surface_integral_of_green_function_partial_xinormal(C(n,:), to.Points(to.ConnectivityList(k,:),:),X,F(n,:));
    end
end

