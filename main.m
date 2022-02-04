% NOTE: FFR -> Flagged For Removal
% NOTE: Replace this with iterative asq???
addpath('asq');

% Makes a traingulation object of the mesh
to = stlread("models/Glider v9.stl");

% For some reason the model is read 10 times too small???
to = triangulation(to.ConnectivityList, to.Points * 10);

% FFR: Plots traingulation object
figure();
trimesh(to);
hold on
axis equal

% Finds each triangle's center point and face normal vectors
CP  = incenter(to);
FN = faceNormal(to);

% FFR: Plotes the triangulation object's normals
quiver3(CP(:,1),CP(:,2),CP(:,3), ...
     FN(:,1),FN(:,2),FN(:,3),0.5,'color','r');

% Records the number of surface triangles
[N,~] = size(to);

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
Ts = [2]; %union(2:10,4:0.2:6); % list of wave periods 
[~,nT] = size(Ts);
thetas = [0,pi()/2,pi()]; % list of incident angles
thetac = ["o","+","*"]; % list of plot modifiers for each theta
[~,ntheta] = size(thetas);

% Defines a list of N triangles so that parfor can nicely slice data
Tri = zeros(3,3,N);
for j = 1:N
    Tri(:,:,j) = to.Points(to.ConnectivityList(j,:),:);
end

M = bodyInertiaMatrix(pm);
C = hydrostaticRestoringMatrix(pm,g);

phis = zeros(nT,N,6);
RAOs = zeros(nT,ntheta,6);

for j = 1:nT
    w = 2*pi()/Ts(j); % sets angular frequency
    K = w^2/g; % sets K(w)


    phis(j,:,:) = velocityPotential(N,CP,FN,FN6,Tri,K);
    [A,B] = addedMassAndDampingMatrices(Tri,phis(j,:,:),FN,p,w);

    for i = 1:ntheta
        [F] = excitingForce(Tri,phis(j,:,:),FN,p,k,g,thetas(i)); 
    
        nu = linsolve(-w^2*(M + A) + 1i*w*B + C,F);
    
        H = nu / a; % response function
        RAOs(j,i,:) = abs(H); % response amplitude operator
    end
end

figure()
hold on
tl = tiledlayout(3,1);

ax1 = nexttile;
title(ax1,'Heave')
ax2 = nexttile;
title(ax2,'Surge')
ax3 = nexttile;
title(ax3,'Sway')

linkaxes([ax1,ax2,ax3],'x');
xlim([0.5 10.5])
ylim([0 inf])
xlabel(tl,'RAO')
ylabel(tl,'Period')
t.TileSpacing = 'compact';

for j = 1:nT
    for i = 1:ntheta
        plot(ax1,Ts(j),RAOs(j,1,3),thetac(i))
        plot(ax2,Ts(j),RAOs(j,1,1),thetac(i))
        plot(ax3,Ts(j),RAOs(j,1,2),thetac(i))
    end
end