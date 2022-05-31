function [VS,CB] = volumeMomentsSubmerged(CL,P,VE)
% Given a 4 element triangulation object abd the volumes of each element
%  returns the submerged volume and center of buoyancy assuming the
%  waterplane is z = 0
    VS = 0;    
    CB = zeros(1,3);

    [N,~] = size(CL);
    for i = 1:N
        tetra = P(CL(i,:),:);

        sgn = sign(tetra(:,3));
        num0 =  sum(sgn(:) == 0); % number of points on z = 0
        sumsgn = sum(sgn);
        if sumsgn == -4 || sumsgn == -3 || (sumsgn == -2 && num0 == 2) ...
        || (sumsgn == -1 && num0 == 3) % all below or at waterline
            tetraCP = sum(tetra)/4;
            tetraV = VE(i);
        elseif sumsgn == -2 || sumsgn == -1 || (sumsgn == 2 && num0 == 0) ...
        || (sumsgn == 1 && num0 == 1) || (sumsgn == 0 && num0 == 2)% three below waterline
            tri = zeros(3,3); % x,y,z cols
            if sumsgn == 0
                sumsgn = -1; % have to remove the positive volume
            end
            [~,loner] = ismember(-sign(sumsgn),sgn); % index of the lone point
            c = 1; % tri index counter
            tetraCP = 0;
            for j = 1:4
                if j ~= loner
                    tri(c,:) = findZeros(tetra(loner,:), tetra(j,:));
                    tetraCP = tetraCP + tetra(j,:) + tri(c,:);
                    c = c + 1;
                end
            end
            tetraCP = tetraCP/6;

            if sumsgn < 0
                % Volume of element tetra - volume of waterplane tetra
                tetraV = VE(i) - volumeTetrahedron([tetra(loner,:);tri]);
            else
                % Volume of underwaterplane tetra
                tetraV = volumeTetrahedron([tetra(loner,:);tri]);
            end
        elseif sumsgn == 0 % two below waterline
            [~,ai1] = ismember(1,sgn); % index of a point above z = 0
            [~,ai2] = ismember(1,sgn(ai1 + 1:4)); % index of a point above
            ai2 = ai2 + ai1;
            [~,bi1] = ismember(-1,sgn); % index of a point below z = 0
            [~,bi2] = ismember(-1,sgn(bi1 + 1:4)); % index of a point below
            bi2 = bi2 + bi1;
            a1 = tetra(ai1,:); % 4 vertices named by location
            a2 = tetra(ai2,:);
            b1 = tetra(bi1,:);
            b2 = tetra(bi2,:);
            p1 = findZeros(a1,b1); % 4 z = 0 intersections named by location
            p2 = findZeros(a2,b1);
            p3 = findZeros(a2,b2);
            p4 = findZeros(a1,b2);
            tetraCP = (p1 + p2 + p3 + p4 + b1 + b2) / 6; % cp of waterline vertices
            tetraV = volumeTetrahedron([b1;p1;p2;p4]) + ...
                     volumeTetrahedron([b1;p2;p3;p4]) + ...
                     volumeTetrahedron([b1;b2;p3;p4]);
        else
            tetraV = 0;
            tetraCP = 0;
        end

        VS = VS + tetraV;
        CB = CB + tetraV*tetraCP;
    end

    CB = CB/VS;    
end