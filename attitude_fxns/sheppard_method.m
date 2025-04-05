function [q] = sheppard_method(C)
%SHEPPARD_METHOD Employs the Sheppard's method to find quaternions from a
%DCM
%   Inputs:
%           C  :  3x3 DCM
%   Outputs:
%           q  :  4x1 quaternion

q1 = sqrt(1/4*(1 + 2*C(1,1) - trace(C)));
q2 = sqrt(1/4*(1 + 2*C(2,2) - trace(C)));
q3 = sqrt(1/4*(1 + 2*C(3,3) - trace(C)));
q4 = sqrt(1/4*(1+trace(C)));

switch max([q1 q2 q3 q4])
    case q1
        q2 = (C(1,2) + C(2,1)) / (4 * q1);
        q3 = (C(3,1) + C(1,3)) / (4*  q1);
        q4 = (C(2,3) - C(3,2)) / (4 * q1);
    case q2
        q1 = (C(1,2) + C(2,1)) / (4 * q2);
        q3 = (C(2,3) + C(3,2)) / (4 * q2);
        q4 = (C(3,1) - C(1,3)) / (4 * q2);
    case q3
        q1 = (C(3,1) + C(1,3)) / (4 * q3);
        q2 = (C(2,3) + C(3,2)) / (4 * q3);
        q4 = (C(1,2) - C(2,1)) / (4 * q3);
    case q4
        q1 = (C(2,3) - C(3,2)) / (4 * q4);
        q2 = (C(3,1) - C(1,3)) / (4 * q4);
        q3 = (C(1,2) - C(2,1)) / (4 * q4);
end
q = [q1; q2; q3 ;q4] / norm([q1; q2; q3 ;q4]);
end