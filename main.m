%function main()
%% Loads And Plots Glider's Mesh
% NOTE: FFR -> Flagged For Removal

% NOTE: Replace this with iterative asq???
addpath('asq');

model = createpde();

% Creates a new pde object for mesh generation
modelpath = "models/Glider v0";
bodypath = " Wingless.stl";
wingpath = " Wings.stl";
importGeometry(model,append(modelpath,bodypath))

% Scales geometry to correct size
scale(model.Geometry, 0.0001);

% FFR: Plots geometry
figure(3);
pdegplot(model,'FaceLabels','on');

% Generates a finite element tetrahedra mesh object of the geometry
generateMesh(model, 'GeometricOrder','linear','Hmin',0.01);

figure(4)
hold on

% Connectivity and point list of waterplane triangles
[CL,P] = waterplaneTriangulation(model,2);

% Plots waterplane triangle points
[n,~] = size(CL);
for i = 1:3
    for j = 1:n
        plot3(P(CL(j,i),1),P(CL(j,i),2),P(CL(j,i),3),"or")
    end
end

% FFR: Plots mesh
pdeplot3D(model.Mesh);
hold off

% Gets waterplane area, center of flotation, and two moments
[S,CF,S11,S22] = waterplaneMoments(CL,P);

% Total Volume
V = volume(model.Mesh);

% Center of Buoyancy
[~,nn] = size(model.Mesh.Nodes);
CB = sum(model.Mesh.Nodes,2)/nn - CF;

% Imports body and wing meshes as triangulation objects
bto = stlread(append(modelpath,bodypath));
wto = stlread(append(modelpath,wingpath));

[nbP,~] = size(bto.Points);
wtoCL = wto.ConnectivityList + [nbP,nbP,nbP];
toCL = bto.ConnectivityList; %[bto.ConnectivityList;wtoCL];
toP = bto.Points*0.0001 - CF.'; %[bto.Points*0.0001;wto.Points*0.001] - CF.';

to = triangulation(toCL, toP);

[N,~] = size(to.ConnectivityList);

% FFR: Plots traingulation object
figure(1);
trimesh(to);
hold on
axis equal

% Finds each triangle's center point and face normal vectors
CP  = incenter(to);
FN = faceNormal(to);

% FFR: Plotes the triangulation object's normals
quiver3(CP(:,1),CP(:,2),CP(:,3), ...
     FN(:,1),FN(:,2),FN(:,3),0.5,'color','r');

% Creates each triangle's 6-dimensional normal vector
FN6 = normal6DOF(CP,FN);

% Plotes the triangulation object's rotational normals
quiver3(CP(:,1),CP(:,2),CP(:,3), ...
     FN6(:,4),FN6(:,5),FN6(:,6),0.5,'color','g');
%% Loads Saved Data

% Defines save files
tstr = append(modelpath,'-ts.txt'); 
raostr = append(modelpath,'-raos.txt'); 
phistr = append(modelpath,'-phis.txt'); 

% Reads from savefiles
rts = readmatrix(tstr); % replace with [] if file doesn't exist
rraos = readmatrix(raostr); % replace with [] if file doesn't exist
rphis = readmatrix(phistr); % replace with [] if file doesn't exist
%% Defines RAO Variables

% Temporarily define response inputs
g = 9.81; 
p = 997; % water density
l = 106.4; % wavelength
k = 2*pi()/l; % wavenumber
a = 1; % wave amplitude
z = @(x) sin(10*pi()/180)*(x - CB(1)) + CB(3) - 0.02;
x = @(a) CB(1) + a;
pm = [[10,(CB.' - [x(-.4),0,z(x(-.4))])];
    [10,(CB.' - [x(-.3),0,z(x(-.3))])];
    [10,(CB.' - [x(-.2),0,z(x(-.2))])];
    [10,(CB.' - [x(-.1),0,z(x(-.1))])];
    [10,(CB.' - [x(.0),0,z(x(.0))])];
    [10,(CB.' - [x(.1),0,z(x(.1))])];
    [10,(CB.' - [x(.2),0,z(x(.2))])];
    [10,(CB.' - [x(.3),0,z(x(.3))])];
    [10,(CB.' - [x(.4),0,z(x(.4))])]];
nts = 1:7; %union(2:10,4:0.2:6); % list of wave periods 
Ts = setdiff(nts,rts); 
[~,nT] = size(Ts);
thetas = [0,pi()/2,pi()]; % list of incident angles
thetac = ["o","+","*"]; % list of plot modifiers for each theta
[~,ntheta] = size(thetas);
%% Solves For The Glider's RAOs

M = bodyInertiaMatrix(pm);
C = hydrostaticRestoringMatrix(pm,g,p,V,CB,S,S11,S22);

% Defines a list of N triangles so that parfor can nicely slice data
Tri = zeros(3,3,N);
for j = 1:N
    Tri(:,:,j) = to.Points(to.ConnectivityList(j,:),:);
end

phis = zeros(6,N,nT);
RAOs = zeros(6,ntheta,nT);
for j = 1:nT
    w = 2*pi()/Ts(j); % sets angular frequency
    K = w^2/g; % sets K(w)

    phis(:,:,j) = velocityPotential(N,CP,FN,FN6,Tri,K);
    [A,B] = addedMassAndDampingMatrices(Tri,phis(:,:,j),FN,p,w);

    for i = 1:ntheta
        [F] = excitingForce(Tri,phis(:,:,j),FN,p,k,g,thetas(i)); 
    
        nu = linsolve(-w^2*(M + A) + 1i*w*B + C,F);
    
        H = nu / a; % response function
        RAOs(:,i,j) = abs(H); % response amplitude operator
    end
end
%% Writes To New Save Files

% Combines saved data with new data
Tcs = [Ts,setdiff(rts,Ts)];
[~,nTc] = size(Tcs);
phics = zeros(6,N,nTc);
phics(:,:,1:nT) = phis(:,:,:);
RAOcs = zeros(6,ntheta,nTc);
RAOcs(:,:,1:nT) = RAOs(:,:,:);
for j = 1:nTc - nT
    phics(:,:,nT+j) = rphis(6*j-5:6*j,:);
    RAOcs(:,:,nT+j) = rraos(6*j-5:6*j,:);
end

% Writes combined data to file
writematrix(Tcs, tstr);
fid = fopen(raostr,'w');
fclose(fid);
fid = fopen(phistr,'w');
fclose(fid);
for j = 1:nTc
    writematrix(RAOcs(:,:,j),raostr,'WriteMode','append')
    writematrix(phics(:,:,j),phistr,'WriteMode','append')
end
%% Plots RAOs

figure(2)
tl = tiledlayout(3,1);

ax1 = nexttile;
title(ax1,'Heave')
hold on
ax2 = nexttile;
title(ax2,'Surge')
hold on
ax3 = nexttile;
title(ax3,'Sway')
hold on

xlim([0.5 10.5])
ylim([0 6])
linkaxes([ax1,ax2,ax3],'xy')
xlabel(tl,'Period')
ylabel(tl,'RAO')
tl.TileSpacing = 'compact';

for j = 1:nTc
    for i = 1:ntheta
        plot(ax1,Tcs(j),RAOcs(3,i,j),thetac(i))
        plot(ax2,Tcs(j),RAOcs(1,i,j),thetac(i))
        plot(ax3,Tcs(j),RAOcs(2,i,j),thetac(i))
    end
end
%end