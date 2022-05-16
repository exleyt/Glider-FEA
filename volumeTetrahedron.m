function [V] = volumeTetrahedron(tetrahedra)
% Returns the volume of a tetrahedra defined as a 4x3 matrix of points
    x1 = tetrahedra(1,1);
    x2 = tetrahedra(2,1);
    x3 = tetrahedra(3,1);
    x4 = tetrahedra(4,1);
    y1 = tetrahedra(1,2);
    y2 = tetrahedra(2,2);
    y3 = tetrahedra(3,2);
    y4 = tetrahedra(4,2);
    z1 = tetrahedra(1,3);
    z2 = tetrahedra(2,3);
    z3 = tetrahedra(3,3);
    z4 = tetrahedra(4,3);

    V = abs((x4 - x1)*((y2 - y1)*(z3 - z1) - (z2 - z1)*(y3 - y1)) + ...
            (y4 - y1)*((z2 - z1)*(x3 - x1) - (x2 - x1)*(z3 - z1)) + ...
            (z4 - z1)*((x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1)))/6;
end

