clc; clear; close all;

script_path = matlab.desktop.editor.getActiveFilename;
run(fullfile(fileparts(mfilename('fullpath')), 'add_to_path.m'));

syms rho theta phi omega_1 omega_2 omega_3 q_sat_1 q_sat_2 q_sat_3 q_sat_4 q_ant_1 q_ant_2 q_ant_3 q_ant_4 e_theta e_phi theta_max phi_max error_weight
syms lambda_1 lambda_2 lambda_3 lambda_4 lambda_5 lambda_6 lambda_7 lambda_8 lambda_9 lambda_10 lambda_11 lambda_12 lambda_13 lambda_14 lambda_15 lambda_16
syms mu_1 mu_2 mu_3 mu_4 mu_5 u_1 u_2 u_3 u_4 u_5
syms I_ant_11 I_ant_12 I_ant_13 I_ant_21 I_ant_22 I_ant_23 I_ant_31 I_ant_32 I_ant_33
%syms I_sat_11 I_sat_12 I_sat_13 I_sat_21 I_sat_22 I_sat_23 I_sat_31 I_sat_32 I_sat_33
syms l m_sat m_ant r_cg_1 r_cg_2 r_cg_3

% form state
w = [omega_1, omega_2, omega_3].';
q_LOS2SAT = [q_sat_1, q_sat_2, q_sat_3, q_sat_4].';
q_SAT2ANT = [q_ant_1, q_ant_2, q_ant_3, q_ant_4].';
err = [e_theta, e_phi].';

x = [rho; theta; phi; w; q_LOS2SAT; q_SAT2ANT; err];

% form control
u = [u_1 u_2 u_3 u_4 u_5].';

% get state length
n = length(x);

% form weight matrices
Q = sym(zeros(n)); Q(n-1:n, n-1:n) = diag(error_weight*sym(ones(2,1)));

% get MOI
I_ant = [I_ant_11 I_ant_12 I_ant_13; I_ant_21 I_ant_22 I_ant_23; I_ant_31 I_ant_32 I_ant_33];
%I_sat = [I_sat_11 I_sat_12 I_sat_13; I_sat_21 I_sat_22 I_sat_23; I_sat_31 I_sat_32 I_sat_33];
I_sat = PMOI_rectangular_prism(m_sat, [l;l;l]); % m is mass of satellite l is length (its a cubesat, probably 50kg later)

DCM = quat2dcm(q_SAT2ANT); % right quaternion?
I_total = I_sat + DCM*parallel_axis_thm([r_cg_1; r_cg_2; r_cg_3], m_ant, I_ant)*DCM.'; % <-- need parellel axis theorem

% form R
R = diag([I_total, I_total, I_total, I_ant, I_ant, I_ant]);

% form constraints
c = [theta_max; -theta_max; phi_max; -phi_max];

% form costate
lambda = [lambda_1 lambda_2 lambda_3 lambda_4 lambda_5 lambda_6 lambda_7 lambda_8 lambda_9 lambda_10 lambda_11 lambda_12 lambda_13 lambda_14 lambda_15 lambda_16].';
mu = [mu_1 mu_2 mu_3 mu_4 mu_5].';

% dynamics !

