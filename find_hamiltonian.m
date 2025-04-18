clc; clear; close all;

script_path = matlab.desktop.editor.getActiveFilename;
run(fullfile(fileparts(mfilename('fullpath')), 'add_to_path.m'));

syms rho theta phi drho dtheta dphi ddrho ddtheta ddphi omega_sat_1 omega_sat_2 omega_sat_3 omega_ant_1 omega_ant_2 omega_ant_3 
syms q_sat_1 q_sat_2 q_sat_3 q_sat_4 q_ant_1 q_ant_2 q_ant_3 q_ant_4 e_theta e_phi theta_max phi_max error_weight real
syms lambda_1 lambda_2 lambda_3 lambda_4 lambda_5 lambda_6 lambda_7 lambda_8 lambda_9 lambda_10 lambda_11 lambda_12 lambda_13 lambda_14 lambda_15 lambda_16 lambda_17 lambda_18 lambda_19 lambda_20 lambda_21 lambda_22
syms mu_1 mu_2 mu_3 mu_4 u_1 u_2 u_3 u_4 u_5
syms I_ant_11 I_ant_12 I_ant_13 I_ant_21 I_ant_22 I_ant_23 I_ant_31 I_ant_32 I_ant_33
%syms I_sat_11 I_sat_12 I_sat_13 I_sat_21 I_sat_22 I_sat_23 I_sat_31 I_sat_32 I_sat_33
syms l m_sat m_ant r_cg_1 r_cg_2 r_cg_3 t

% form state
w_sat = [omega_sat_1, omega_sat_2, omega_sat_3].';
w_ant = [omega_ant_1, omega_ant_2, omega_ant_3].';
q_LOS2SAT = [q_sat_1, q_sat_2, q_sat_3, q_sat_4].';
q_SAT2ANT = [q_ant_1, q_ant_2, q_ant_3, q_ant_4].';
err = [e_phi, e_theta].';


%% ORDER IS RHO, PHI, THETA %%
x = [rho; phi; theta; drho; dphi; dtheta; w_sat; w_ant; q_LOS2SAT; q_SAT2ANT; err]; 

% form control
u = [u_1 u_2 u_3 u_4 u_5].';

% get state length
n = length(x);

% form weight matrices
Q = sym(zeros(n)); Q(n-1:n, n-1:n) = diag(error_weight*sym(ones(2,1)));

% get MOI
I_ant = [I_ant_11 I_ant_12 I_ant_13; I_ant_21 I_ant_22 I_ant_23; I_ant_31 I_ant_32 I_ant_33]; % <--- needs to be a function based on geometry
I_sat = PMOI_rectangular_prism(m_sat, [l;l;l]); % m is mass of satellite l is length (its a cubesat, probably 50kg later)

DCM = quat2dcm(q_SAT2ANT); % right quaternion?
I_total = I_sat + DCM*parallel_axis_thm([r_cg_1; r_cg_2; r_cg_3], m_ant, I_ant)*DCM.'; % <-- need parellel axis theorem

% form R
%R = diag([I_total, I_total, I_total, I_ant, I_ant, I_ant]);
R = diag([trace(I_total), trace(I_total), trace(I_total), trace(I_ant), trace(I_ant)]);

% form constraints
c = [theta_max; -theta_max; phi_max; -phi_max];

% form costate
lambda = [lambda_1 lambda_2 lambda_3 lambda_4 lambda_5 lambda_6 lambda_7 lambda_8 lambda_9 lambda_10 lambda_11 lambda_12 lambda_13 lambda_14 lambda_15 lambda_16 lambda_17 lambda_18 lambda_19 lambda_20 lambda_21 lambda_22].';
mu = [mu_1 mu_2 mu_3 mu_4].';

%% dynamics !

% get vel and acc
v = [drho; rho*dphi; rho*dtheta*sin(phi)];
a = [
    (ddrho - rho*dphi^2 - rho*dtheta^2 * sin(phi)^2);
    rho*ddphi + 2*drho*dphi - rho*dtheta^2*sin(phi)*cos(phi);
    rho*ddtheta*sin(phi) + 2*drho*dtheta*sin(phi) + 2*rho*dphi*dtheta*cos(phi)
];

% get wumbo rate of satellite body wrt inertial (note order of u_s to axis # !)
% recall dw wants inertia tensor diagonal soo
dw_sat = [
    (-(I_total(3,3) - I_total(2,2)) * omega_sat_2*omega_sat_3 + u_1) / I_total(1,1);
    (-(I_total(1,1) - I_total(3,3)) * omega_sat_3*omega_sat_1 + u_2) / I_total(2,2);
    (-(I_total(2,2) - I_total(1,1)) * omega_sat_1*omega_sat_2 + u_3) / I_total(3,3);
];

% get quanternion rate LOS2SAT
dq_LOS2SAT = quat_kde(q_LOS2SAT, w_sat);

% get antenna body w wrt sat body (note a_1 a_2 is roll. let 
dw_ant = [
    (-(I_ant(3,3) - I_ant(2,2)) * omega_ant_2*omega_ant_3) / I_ant(1,1);
    (-(I_ant(1,1) - I_ant(3,3)) * omega_ant_3*omega_ant_1 + u_4) / I_ant(2,2);
    (-(I_ant(2,2) - I_ant(1,1)) * omega_ant_1*omega_ant_2 + u_5) / I_ant(3,3);
];

% quaternion rate SAT2ANT
dq_SAT2ANT = quat_kde(q_SAT2ANT, w_ant); 

% get DCM LOS2ANT 
DCM_LOS2SAT = quat2dcm(q_LOS2SAT);
DCM_SAT2ANT = quat2dcm(q_SAT2ANT);
DCM_LOS2ANT = DCM_LOS2SAT*DCM_SAT2ANT;

% get w_LOS2ANT
omega_LOS2ANT = DCM_SAT2ANT*w_sat + w_ant;


% from my diagram, 2 is roll. 1 and 3 are pitch or yaw, u choose chande
% let 1 be pitch=theta and 3=yaw=phi --> use 31 E-A
% DCM_31 = R3(e_phi) * R1(e_theta); % <-- jut for seeing relationships
e_phi = acos(DCM_LOS2ANT(1,1));
e_theta = acos(DCM_LOS2ANT(3,3));


% finally form the state vector
dx = [
  v;
  a;
  dw_sat;
  dw_ant;
  dq_LOS2SAT;
  dq_SAT2ANT;
  e_phi;
  e_theta
];


% hamiltonian
H = x.'*Q*x + u.'*R*u + lambda.'*dx + mu.'*c;
latex(simplify(H))





% pitch is theta, yaw is phi (this is needed for later, but cant use in
% syms-land)
%[e_theta, e_phi, ~] = dcm2angle(DCM_LOS2ANT, 'ZYX');

