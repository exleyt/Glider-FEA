function [p] = findZeros(p1,p2)
% Returns the point where the line crosses z = 0 (errors if the line is
%  parallel to z = 0)
    p = zeros(1,3);
    t = -p1(3)/(p2(3) - p1(3)); % z = 0
    p(1:2) = (p2(1:2) - p1(1:2))*t + p1(1:2);
end

