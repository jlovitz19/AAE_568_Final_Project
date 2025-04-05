function L = gravity_induced_torque(mu, R, I)
%GRAVITY_INDUCED_TORQUE returns the gravity gradient induced torque up to
% order(epsilon^2) where epsilon = r/R_c
%   inputs:
%       mu      :  orbital parameter
%       R       :  radius [3x1]
%       I       :  intertia tensor
%   outputs:
%       L_g      : gravity gradient induced torque

L = 3*mu / norm(R)^5 * (cross_mat(R) * I) * R;
end

