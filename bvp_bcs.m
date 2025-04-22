function bcs = bvp_bcs(yi,yf,param)
% Extract parameters
sphere0 = param.sphere0;
q_LOS2SAT0 = param.q_LOS2SAT0;
q_SAT2ANT0 = param.q_SAT2ANT0;
e_phi0 = param.e_phi0;
e_theta0 = param.e_theta0;
x_np_1 = param.x_np_1;
lambda0 = param.lambda0;

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

% Lambda
bcs_lambda = [
    yf(33) - lambda0(7);
    yf(34) - lambda0(8);
    yf(35) - lambda0(9);
    yf(36) - lambda0(10);
    yf(37) - lambda0(11);
    yf(38) - lambda0(12);
];

% Initial quaternions
bcs_quat = [
    % Initial
    yi(13) - q_LOS2SAT0(1);
    yi(14) - q_LOS2SAT0(2);
    yi(15) - q_LOS2SAT0(3);
    yi(16) - q_LOS2SAT0(4);
    yi(17) - q_SAT2ANT0(1);
    yi(18) - q_SAT2ANT0(2);
    yi(19) - q_SAT2ANT0(3);
    yi(20) - q_SAT2ANT0(4);
    % Initial
    yf(13) - q_LOS2SAT0(1);
    yf(14) - q_LOS2SAT0(2);
    yf(15) - q_LOS2SAT0(3);
    yf(16) - q_LOS2SAT0(4);
    yf(17) - q_SAT2ANT0(1);
    yf(18) - q_SAT2ANT0(2);
    yf(19) - q_SAT2ANT0(3);
    yf(20) - q_SAT2ANT0(4);
];

% Error
bcs_e = [
    % Initial
    yi(21) - e_phi0;
    yi(22) - e_theta0;
    % Final
    yf(21) - e_phi0;
    yf(22) - e_theta0;
];

% Extra states
bcs_x_np1 = [
    % Initial
    yi(23) - x_np_1(1);
    yi(24) - x_np_1(2);
    yi(25) - x_np_1(3);
    yi(26) - x_np_1(4);
    % Final
    yf(23) - x_np_1(1);
    yf(24) - x_np_1(2);
    yf(25) - x_np_1(3);
    yf(26) - x_np_1(4);
];

bcs = [bcs_sphere; bcs_lambda; bcs_quat; bcs_e; bcs_x_np1];
end