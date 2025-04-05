function C = mrp2DCM(sigma)
%MRP2DCM Given a vector of MRPs returns a DCM
%   Detailed explanation goes here

sig_cross = cross_mat(sigma);
sig_dot = dot(sigma, sigma);

C = eye(3) + (8*sig_cross*sig_cross - 4*(1 - sig_dot)*sig_cross) / (1 + sig_dot)^2;
end

