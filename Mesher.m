model = createpde();
importGeometry(model,"test_models\Glider.stl");
generateMesh(model, 'GeometricOrder','linear','Hmin',0.01,'Hmax',0.04);
CL = model.Mesh.Elements.';
P0 = model.Mesh.Nodes.'; 

p = 1026;

[~,N] = size(model.Mesh.Elements);
B = p*volume(model.Mesh); % Total volume
VE = zeros(1,N); % List of volumes for each element
for j = 1:N
    VE(j) = volume(model.Mesh,j);
end

pointMasses = readmatrix("test_models\Glider Mass.txt");
[M,CM] = massMoments(pointMasses);

VT = M / p; % target volume

%should fail if target > total
pitch = 0;
EP = 1E-6; % 
E = -1;
%[~,CB] = volumeMomentsSubmerged(CL,P,VE);
[PN,~,CBN] = moveDepth(CL,P0,VE,VT,EP); 
dZ = PN(1,3) - P0(1,3);
CMN = CM + [0,0,dZ];
dXN = CBN(1) - CMN(1);

while abs(dXN) > EP
    dPY = sign(dXN)*10^E*pi/2; % delta pitch
    pitch = pitch + dPY; % total pitch
    roty = [cos(dPY),     0,      sin(dPY);
            0,              1,      0;
            -sin(dPY),    0,      cos(dPY)];
    P = (PN - CMN)*roty + CMN; % Center points on origin about CM and rotate
    [PN,~,CBN] = moveDepth(CL,P,VE,VT,EP);  % Find correct depth V*p = M
    dZ = PN(1,3) - P(1,3); % record change in depth
    CMN = CMN + [0,0,dZ]; % move CM to correct dpeth
    dX = dXN;
    dXN = CBN(1) - CMN(1); 
    if sign(dXN) ~= sign(dX)
        E = E - 1;
    end
end

to = triangulation(CL,PN);
[T,P1] = freeBoundary(to);
mesh = triangulation(T,P1);
figure(1) 

trimesh(mesh);
axis equal