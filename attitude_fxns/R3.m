function [mat] = R3(angle)
%R3 Returns the direction cosine matrix (DCM) about axis 3.
%   C = R3(angle) calculates the body-frame DCM abo
mat = [cos(angle) sin(angle) 0;
    -sin(angle) cos(angle) 0;
    0 0 1];
end