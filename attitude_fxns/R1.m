function [mat] = R1(angle)
%R1 Summary of this function goes here
%   Detailed explanation goes here

mat = [1 0 0;
    0 cos(angle) sin(angle)
    0 -sin(angle) cos(angle)];
end