function [dw] = euler_eom(w, I, L)
%EULER_EOM calculates omega_dot using euler's rotational eom for a rigid
%body under torque
%   Detailed explanation goes here
dw = I\(-cross_mat(w) * I*w + L);
end
