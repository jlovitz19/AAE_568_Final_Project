function [v_cross] = cross_mat(v)
%CROSS_MAT returns the cross product matrix of a 3x1 or 1x3 vector
%   Detailed explanation goes here
v_cross = [0 -v(3) v(2);
    v(3) 0 -v(1);
    -v(2) v(1) 0];
end

