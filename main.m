% Creates a new pde object
model = createpde();

% Generates a discrete geometry object from stl file
importGeometry(model,"template-model\files\CORNER.STL");

% Generates a finite element tetrahedra mesh object of the geometry
generateMesh(model, 'GeometricOrder','linear','Hmin',1);

% Plots mesh 
%pdeplot3D(model);

% Makes a traingulation object of mesh
mto = triangulation(model.Mesh.Elements.', model.Mesh.Nodes.');

% Makes a traingle connectivity and point list of surface triangles
[T, P] = freeBoundary(mto);

% Makes a traingulation object of mesh surface
to = triangulation(T, P);

% Plots traingulation object
trimesh(t);

% Sets N equal to the number of surface triangles
[N,~] = size(T);

% Declares the pde as a system of equations N equations
model.PDESystemSize = N;