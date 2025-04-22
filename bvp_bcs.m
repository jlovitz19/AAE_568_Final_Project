%% ALERT WE ARE STILL 2 CONSTRAINTS SHORT

function bcs = bvp_bcs(yi,yf,param)
% Extract parameters
sphere0 = param.sphere0;
q_LOS2SAT0 = param.q_LOS2SAT0;
q_SAT2ANT0 = param.q_SAT2ANT0;
e_phi0 = param.e_phi0;
e_theta0 = param.e_theta0;
x_np1_0 = param.x_np1_0;

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

bcs_w = [
    % Satellite initial = final
    yi(7)-yf(7);
    yi(8)-yf(8);
    yi(9)-yf(9);
    % Antenna initial = final
    yi(10)-yf(10);
    yi(11)-yf(11);
    yi(12)-yf(12);
];

% Lambda
bcs_lambda = [
    % Lambdas corresponding to angular velocities time
    yf(33);
    yf(34);
    yf(35);
    yf(36);
    yf(37);
    yf(38);
    % Lambdas corresponding to extra states
    yf(49);
    yf(50);
    yf(51);
    yf(52);

    % lambda corresponding to error
    yi(47);
    yi(48);
];

% Initial quaternions
bcs_quat = [
    % Initial --> MIGHT BE OVERCONSRAINING QUATERNIONS
    yi(13) - q_LOS2SAT0(1);
    yi(14) - q_LOS2SAT0(2);
    yi(15) - q_LOS2SAT0(3);
    %yi(16) - q_LOS2SAT0(4);
    yi(17) - q_SAT2ANT0(1);
    yi(18) - q_SAT2ANT0(2);
    yi(19) - q_SAT2ANT0(3);
    %yi(20) - q_SAT2ANT0(4);

    % Normalization constraint
    yi(13)^2+yi(14)^2+yi(15)^2+yi(16)^2 - 1;
    yi(17)^2+yi(18)^2+yi(19)^2+yi(20)^2 - 1;
    yf(13)^2+yf(14)^2+yf(15)^2+yf(16)^2 - 1;
    yf(17)^2+yf(18)^2+yf(19)^2+yf(20)^2 - 1;
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
    yi(23) - x_np1_0(1);
    yi(24) - x_np1_0(2);
    yi(25) - x_np1_0(3);
    yi(26) - x_np1_0(4);
    % Final
    yf(23) - x_np1_0(1);
    yf(24) - x_np1_0(2);
    yf(25) - x_np1_0(3);
    yf(26) - x_np1_0(4);
];

bcs = [bcs_sphere; bcs_w; bcs_lambda; bcs_quat; bcs_e; bcs_x_np1];
end