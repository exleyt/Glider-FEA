% NOTE: FFR -> Flagged For Removal
% NOTE: Replace this with iterative asq???
addpath('asq');

% Makes a traingulation object of the mesh
to = stlread("models/Glider v8.stl");

% FFR: Plots traingulation object
figure(10);
trimesh(to);
hold on
axis equal

% Records the number of surface triangles
[N,~] = size(to);

% Finds each triangle's center point and face normal vectors
CP  = incenter(to);
FN = faceNormal(to);

% FFR: Plotes the triangulation object's normals
quiver3(CP(:,1),CP(:,2),CP(:,3), ...
     FN(:,1),FN(:,2),FN(:,3),0.5,'color','r');

% Creates each triangle's 6-dimensional normal vector
FN6 = zeros(N,3);
FN6(:,1:3) = FN(:,:);
FN6(:,4:6) = cross(CP(:,:),FN(:,:));

% Plotes the triangulation object's rotational normals
quiver3(CP(:,1),CP(:,2),CP(:,3), ...
     FN6(:,4),FN6(:,5),FN6(:,6),0.5,'color','g');

% Temporarily define response inputs
g = 9.81; 
p = 997; % water density
l = 76; % wavelength
k = 2*pi()/l; % wavenumber
a = 4; % wave amplitude
pm = [9,0,0,-0.3
    18,0,0.34,-0.24004882;
    18,0,0.68,-0.18009764;
    18,0,1.2,-0.12014646;
    18,0,1.4,-0.531422;
    9,0,1.7,0]; % point masses along body length (x)

% Defines a list of N triangles so that parfor can nicely slice data
Tri = zeros(3,3,N);
for j = 1:N
    Tri(:,:,j) = to.Points(to.ConnectivityList(j,:),:);
end

M = bodyInertiaMatrix(pm);
C = hydrostaticRestoringMatrix(pm,g);

for T = 1:10
    w = 2*pi()/T;
    K = w^2/g;

    phi = velocityPotential(N,CP,FN,FN6,Tri,K);
    [A,B] = addedMassAndDampingMatrices(Tri,phi,FN,p,w);

    
    angle = [0,pi()/2,pi()]; % incident wave angle
    shape = ["o","+","*"]; % shape of point at for each angle

    % Loops through different incident wave angles
    for j = 1:3
        [F] = excitingForce(Tri,phi,FN,p,k,g,angle(j));
        
        nu = linsolve(-w^2*(M + A) + 1i*w*B + C,F);
        H = nu / a; % response function
        RAO = abs(H); % response amplitude operator
   
        % Plots surge RAO on fig1
        figure(1)
        hold on
        plot(T, RAO(1), shape(j));
        
        % Plots sway RAO on fig2
        figure(2)
        hold on
        plot(T, RAO(2), shape(j));

        % Plots heave RAO on fig3
        figure(3)
        hold on
        plot(T, RAO(3), shape(j));
    end
end