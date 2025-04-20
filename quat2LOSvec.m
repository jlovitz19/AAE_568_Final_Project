function los_vec = quat2LOSvec(LOS, q_LOS2SAT, q_SAT2ANT)
%QUAT2LOSVEC Summary of this function goes here
%   qiven quaternions (or state, could change fxn depending on use-case),
%   find and return the unit vector for pointing and evaluation reasons

R_LOS2SAT = quat2dcm(q_LOS2SAT);
R_SAT2ANT = quat2dcm(q_SAT2ANT);




end

