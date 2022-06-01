function [rot] = roty(angle)
% Returns a 3D rotation matrix about the y axis by angle radians
    rot = [cos(angle),     0,      sin(angle);
           0,              1,      0;
           -sin(angle),    0,      cos(angle)];
end

