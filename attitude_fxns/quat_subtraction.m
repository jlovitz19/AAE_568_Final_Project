function [diff] = quat_subtraction(A,B)
%QUAT_SUBTRACTION computes quaternion A - B
%   Detailed explanation goes here

mat = [B(4) -B(3) B(2) B(1);
    B(3) B(4) -B(1) B(2);
    -B(2) B(1) B(4) B(3);
    -B(1) -B(2) -B(3) B(4)];

diff = mat.' * A;
end

