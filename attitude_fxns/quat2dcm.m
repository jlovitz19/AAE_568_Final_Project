function [C] = quat2dcm(q)
%QUAT2DCM converts quaternions to a DCM
%   inputs: q = quaternions as a 4x1 or 1x4 vector
%   outputs: C = 3x3 DCM
q_cross = cross_mat(q(1:3));
C = eye(3) - 2*q(4) * q_cross + 2*q_cross*q_cross;
end