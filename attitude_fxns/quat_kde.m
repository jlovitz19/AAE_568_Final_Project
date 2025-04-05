function [dq] = quat_kde(q, w)
%QUAT_KDE returns q' given w and q
%   inputs:
%           q  :  set of quaternions
%           w  :  angular velocities
%   outputs:
%           dq :  set of quaternion time derivatives

% q' = 1/2 [\omega matrix] * q  (could also use q matrix *w]

% note top 3x3 of w_mat is cross product matrix transposed
w_mat = [cross_mat(w).', [w(1); w(2); w(3)]; [-w(1) -w(2) -w(3) 0]];

% write q explicity in case a 1x3 q was given
dq = 1/2 * w_mat * [q(1); q(2); q(3); q(4)]; 
end