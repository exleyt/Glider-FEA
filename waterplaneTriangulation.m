function [CL,P2] = waterplaneTriangulation(model,face)
% Find the triangles that make up the waterplane
%
% Takes a meshed model and the face that makes up its waterplane and
%  returns the triangles that belong to that face as a connectivity list 
%  CL(N,3) and list of points P(n,3) that can be used to create a 
%  triangulation object 
    % Finds the nodes of the waterplane
    L1 = findNodes(model.Mesh,'Region','Face',face); % List of node indices
    [~,nL1] = size(L1);
    P = zeros(nL1,3); % List of node points
    for j = 1:nL1
        P(j,:) = model.Mesh.Nodes(:,L1(j));
    end

    % Creates triangulation of mesh of tetrahedra
     mto = triangulation(model.Mesh.Elements.', model.Mesh.Nodes.');
    
    % Makes a traingle connectivity and point list of surface triangles
    [CL2,P2] = freeBoundary(mto);
    
    % Creates a list of waterplane node locations in P2
    [nP2,~] = size(P2);
    LP2 = zeros(nL1, 1);
    for i = 1:nL1
        for j = 1:nP2
            if P(i,:) == P2(j,:)
                LP2(i) = j;
                break;
            end
        end
    end
    
    % Gets length of CL
    [nL2,~] = size(CL2);
    N = 0;
    for j = 1:nL2
        if sum(ismember(CL2(j,:),LP2)) == 3
            N = N + 1;
        end
    end

    CL = zeros(N,3);
    i = 1;
    for j = 1:nL2
        if sum(ismember(CL2(j,:),LP2)) == 3
            CL(i,:) = CL2(j,:);
            i = i + 1;
        end
    end
end

