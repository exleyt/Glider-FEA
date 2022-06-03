function [PN,CB] = moveDepth(CL,P,VE,VT,EP)
% Moves the object mesh (CL,P) vertically until the volume below z = 0
%  is equal to VT. VE is an ordered list of volumes of each element in CL 
%  and EP is the target accuracy/epsilon.
% Returns the new set of points PN and center of buoyancy CB
    E = 0;
    [VS,CB] = volumeMomentsSubmerged(CL,P,VE);

    % New/current loop values
    dVN = VS - VT; % negative delta volume
    PN = P;

    while abs(dVN) > EP
        PN = PN + sign(dVN)*[0,0,1*10^E]; % translate by 1*10^E
        [VSN,CB] = volumeMomentsSubmerged(CL,PN,VE);
        if VSN == 0
            VSN = VT - 1; % will make dVN negative to move object down
        end
        dV = dVN;
        dVN = VSN - VT;
        if sign(dVN) ~= sign(dV)
            E = E - 1;
        end
    end
end

