function bcs = bvp_bcs(yi,yf,param)
% Extract parameters
sphere0 = param.sphere0;
w0 = param.w0;
q_LOS2SAT0 = param.q_LOS2SAT0;
q_SAT2ANT0 = param.q_SAT2ANT0;
e_phi0 = param.e_phi0;
e_theta0 = param.e_theta0;

% Position and velocity in spherical coordinates
bcs_sphere = [
    % Initial
    yi(1) - sphere0(1);
    yi(2) - sphere0(2);
    yi(3) - sphere0(3);
    yi(4) - sphere0(4);
    yi(5) - sphere0(5);
    yi(6) - sphere0(6);
    % Final
    yf(1) - sphere0(1);
    yf(2) - sphere0(2);
    yf(3) - sphere0(3);
    yf(4) - sphere0(4);
    yf(5) - sphere0(5);
    yf(6) - sphere0(6);
];

% Angular velocities
bcs_w = [
    yi(7) - w0(1);
    yi(8) - w0(2);
    yi(9) - w0(3);
    yi(10) - w0(4);
    yi(11) - w0(5);
    yi(12) - w0(6);
];

% Initial quaternions
bcs_quat = [
    yi(13) - q_LOS2SAT0(1);
    yi(14) - q_LOS2SAT0(2);
    yi(15) - q_LOS2SAT0(3);
    yi(16) - q_LOS2SAT0(4);
    yi(17) - q_SAT2ANT0(1);
    yi(18) - q_SAT2ANT0(2);
    yi(19) - q_SAT2ANT0(3);
    yi(20) - q_SAT2ANT0(4);
];

% Initial error
bcs_e = [
    yi(21) - e_phi0
    yi(22) - e_theta0;
];


%{
bcs = [
    % Spherical coords
    yi(1)-sphere0(1), yi(2)-sphere0(2),... 
    yi(3)-sphere0(3), yi(4)-sphere0(4), yi(5)-sphere0(5), yi(6)-sphere0(6),...
    % Angular velocities
    yi(7)-w0(1), yi(8)-w0(2), yi(9)-wo(3),... 
    yi(10)-w0(4), yi(11)-w0(5), yi(12)-w0(6),... 

    yi(13), yi(14), yi(15), yi(16), yi(17), yi(18), yi(19), yi(20),...
    yi(21), yi(22), 

    yf(1)-yi(1), yf(2)-yi(2), yf(3)-yi(3),...
    yf(4)-yi(4), yf(5)-yi(5), yf(6)-yi(6),...

    yf(7), yf(8), yf(9), yf(10), yf(11), yf(12), yf(13),...
    yf(14), yf(15), yf(16), yf(17), yf(18), yf(19), yf(20), yf(21),...
    yf(22), yf(23), yi(23), yf(24), yi(24), yf(25), yi(25),...
    yf(26), yi(26)
];
%}
bcs = [bcs_sphere; bcs_w; bcs_quat; bcs_e
end