function [mat] = R2(angle)
%R2 Summary of this function goes here
%   Detailed explanation goes here
mat = [cos(angle) 0 -sin(angle);
    0 1 0;
    sin(angle) 0 cos(angle)];
end