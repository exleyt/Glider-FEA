function [V,CB] = volumeMoments(mesh)
% Returns the volume and center of buoyancy of a mesh
%
% Given a mesh returns: 
% V the submerged volume
% CB the volume's center
    V = volume(mesh);
    CB = zeros(3,1);

    [~,N] = size(mesh.Elements);

    for j = 1:N
        e = mesh.Elements(:,j);
        tetra = mesh.Nodes(:,e);
        tetraCP = sum(tetra,2)/4;
        tetraV = volume(mesh,e);

        CB = CB + tetraV*tetraCP;
    end

    CB = CB/V;
end