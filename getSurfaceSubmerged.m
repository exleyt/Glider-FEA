function [TO] = getSurfaceSubmerged(CL,P,EP)
% Given a 3 element triangulation object returns the submerged surface
%  assuming the waterplane is z = 0
    [nP,~] = size(P);
    PE = zeros(2*nP,3);
    PE(1:nP,:) = P; % array for extra points
    pi = nP;
    [nCL,~] = size(CL);
    CLE = zeros(3*nCL,3); % array for extra connections
    cli = 0;

    for i = 1:nCL
        tri = P(CL(i,:),:);
        sgn = sign(tri(:,3));
        num0 = sum(ismember(0,tri(:,3))); % number of points on z = 0
        sumsgn = sum(sgn);
        if sumsgn == -3 || (sumsgn == -2 && num0 == 1) ...
        || (sumsgn == -1 && num0 == 2) % all below or at waterline
            cli = cli + 1;
            CLE(cli,:) = CL(i,:);
        elseif sumsgn == -1 || (sumsgn == 1 && num0 == 0) % 2 above/below waterline
            [~,ai] = ismember(-sign(sumsgn),sgn); % index of the lone point
            bi = mod(ai,3) + 1;
            ci = mod(bi,3) + 1;
            a = tri(ai,:);
            b = tri(bi,:);
            c = tri(ci,:);
            d = findZeros(a,b);
            e = findZeros(a,c);
            di = -1; % index of d in PE, not tri/CL like a,b,c
            ei = -1; % index of e in PE, not tri/CL like a,b,c
            for j = nP + 1:pi % loops over all added PEs
                if norm(d - PE(j,:)) < EP
                    di = j;
                elseif norm(e - PE(j,:)) < EP
                    ei = j;
                end
            end
            if di == -1
                pi = pi + 1;
                PE(pi,:) = d;
                di = pi;
            end
            if ei == -1
                pi = pi + 1;
                PE(pi,:) = e;
                ei = pi;
            end

            if sumsgn == -1 
                cli = cli + 1;
                CLE(cli,:) = [CL(i,bi),CL(i,ci),di];
                cli = cli + 1;
                CLE(cli,:) = [CL(i,ci),ei,di]; % order matters for normals
            else
                cli = cli + 1;
                CLE(cli,:) = [CL(i,ai),di,ei];
            end           
        elseif sumsgn == 0 && num0 == 1
            [~,ai] = ismember(0,sgn); % index of the z = 0 point
            bi = mod(ai,3) + 1;
            ci = mod(bi,3) + 1;
            b = tri(bi,:);
            c = tri(ci,:);
            d = findZeros(b,c);
            di = -1; % index of d in PE, not tri/CL like a,b,c
            for j = nP + 1:pi % loops over all added PEs
                if norm(d - PE(j,:)) < EP
                    di = j;
                end
            end
            if di == -1
                pi = pi + 1;
                PE(pi,:) = d;
                di = pi;
            end
            if sgn(bi) == -1
                cli = cli + 1;
                CLE(cli,:) = [CL(i,ai),CL(i,bi),di];
            else
                cli = cli + 1;
                CLE(cli,:) = [CL(i,ai),di,CL(i,ci)];
            end
        end
    end

    CLFN = CLE(1:cli,:); % array of every triangle with maps to PE    
    PFNM = zeros(pi,1); % array that maps used points in PFN
    for i = 1:cli
        PFNM(CLFN(i,:)) = [1,1,1];
    end

    PN = zeros(sum(PFNM),3); % array of all used points
    PEM = zeros(pi,1); % maps each point in PE to a point in PN
    c = 0;
    for i = 1:pi
        if PFNM(i) == 1
            c = c + 1;
            PN(c,:) = PE(i,:);
            PEM(i) = c;
        end
    end

    CLN = zeros(cli,3);
    for i = 1:cli
        CLN(i,:) = PEM(CLFN(i,:));
    end

    try
        TO = triangulation(CLN,PN);
    catch
        TO = [];
    end
end 