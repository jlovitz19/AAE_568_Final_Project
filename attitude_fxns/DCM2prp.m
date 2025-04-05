function [lambda,theta] = DCM2prp(C)
%EULER_ANGLES2PRP converts a DCM to principal rotation parameters lambda
%and theta
%   Detailed explanation goes here

theta = acos( (trace(C) - 1) / 2 );
lambda = 1/(2*sin(theta)) * [C(2,3) - C(3,2); C(3,1) - C(1,3); C(1,2) - C(2,1)];
end