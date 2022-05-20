function [TO] = getSurfaceSubmerged(CL,P,EP)
% Given a 3 element triangulation object returns the submerged surface
%  assuming the waterplane is z = 0
    [nP,~] = size(P);
    PE = zeros(2*nP,3);
    PE(1:nP,:) = P;
    pi = nP + 1;
    [nCL,~] = size(CL);
    CLE = zeros(3*nCL,3);
    cli = 1;

    for i = 1:nCL
        tri = P(CL(i,:),:);
        sgn = sign(tri(:,3));
        num0 = sum(ismember(0,tri(:,3))); % number of points on z = 0
        sumsgn = sum(sgn);
        if sumsgn == -3 || (sumsgn == -2 && num0 == 1) ...
        || (sumsgn == -1 && num0 == 2) % all below or at waterline
            CLE(cli,:) = CL(i,:);
            cli = cli + 1;
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
            for j = nP + 1:pi - 1 % loops over all added PEs
                if norm(d - PE(j,:)) < EP
                    di = j;
                elseif norm(e - PE(j,:)) < EP
                    ei = j;
                end
            end
            if di == -1
                PE(pi,:) = d;
                di = pi;
                pi = pi + 1;
            end
            if ei == -1
                PE(pi,:) = e;
                ei = pi;
                pi = pi + 1;
            end

            if sumsgn == -1 
                CLE(cli,:) = [CL(i,bi),CL(i,ci),di];
                cli = cli + 1;
                CLE(cli,:) = [CL(i,ci),ei,di]; % order matters for normals
                cli = cli + 1;
            else            
                CLE(cli,:) = [CL(i,ai),di,ei];
                cli = cli + 1;
            end           
        elseif sumsgn == 0 && num0 ~= 3
            [~,ai] = ismember(0,sgn); % index of the z = 0 point
            bi = mod(ai,3) + 1;
            ci = mod(bi,3) + 1;
            b = tri(bi,:);
            c = tri(ci,:);
            d = findZeros(b,c);
            di = -1; % index of d in PE, not tri/CL like a,b,c
            for j = nP + 1:pi - 1 % loops over all added PEs
                if norm(d - PE(j,:)) < EP
                    di = j;
                end
            end
            if di == -1
                PE(pi,:) = d;
                di = pi;
                pi = pi + 1;
            end
            CLE(cli,:) = [CL(i,ai),CL(i,bi),di];
            cli = cli + 1;
            CLE(cli,:) = [CL(i,ai),di,CL(i,ci)];
            cli = cli + 1;

        end
    end
    CLN = CLE(1:cli-1,:);
    PN = PE(1:pi-1,:);
    try
        TO = triangulation(CLN,PN);
    catch
        TO = [];
    end
end 