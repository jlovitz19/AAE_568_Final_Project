clc; clear; close all;

script_path = matlab.desktop.editor.getActiveFilename;
run(fullfile(fileparts(mfilename('fullpath')), 'add_to_path.m'));

syms rho theta phi drho dtheta dphi ddrho ddtheta ddphi omega_sat_1 omega_sat_2 omega_sat_3 omega_ant_1 omega_ant_2 omega_ant_3 
syms q_sat_1 q_sat_2 q_sat_3 q_sat_4 q_ant_1 q_ant_2 q_ant_3 q_ant_4 e_theta e_phi error_weight real
syms lambda_1 lambda_2 lambda_3 lambda_4 lambda_5 lambda_6 lambda_7 lambda_8 lambda_9 lambda_10 lambda_11 lambda_12 lambda_13 lambda_14 lambda_15 lambda_16 lambda_17 lambda_18 lambda_19 lambda_20 lambda_21 lambda_22
syms mu_1 mu_2 mu_3 mu_4 u_1 u_2 u_3 u_4 u_5 
assume(q_sat_1^2 + q_sat_2^2 + q_sat_3^2 + q_sat_4^2 == 1)
assume(q_ant_1^2 + q_ant_2^2 + q_ant_3^2 + q_ant_4^2 == 1)

l = .5;
m_sat = 50;
r_max = 0.25;
d_max = 0.05;
t_max = 0.01;
m_ant = 11.6;

% form state
w_sat = [omega_sat_1, omega_sat_2, omega_sat_3].';
w_ant = [omega_ant_1, omega_ant_2, omega_ant_3].';
q_LOS2SAT = [q_sat_1, q_sat_2, q_sat_3, q_sat_4].';
q_SAT2ANT = [q_ant_1, q_ant_2, q_ant_3, q_ant_4].';
err = [e_phi, e_theta].';


%% ORDER IS RHO, PHI, THETA %%

% erm are dphi dro and dtheta supposed ?
x = [rho; phi; theta; drho; dphi; dtheta; w_sat; w_ant; q_LOS2SAT; q_SAT2ANT; err]; 

% form control
u = [u_1 u_2 u_3 u_4 u_5].';

% get state length
n = length(x);

% form weight matrices
Q = sym(zeros(n)); Q(n-1:n, n-1:n) = diag(error_weight*sym(ones(2,1)));

% get MOI
[I_ant, r_cg_ant] = ant_inertia(r_max, d_max, t_max, m_ant);
I_sat = PMOI_rectangular_prism(m_sat, [l;l;l]); % m is mass of satellite l is length (its a cubesat, probably 50kg later)

DCM = quat2dcm(q_SAT2ANT); % right quaternion?
r_cg = r_cg_ant  + ones(3,1) * l / 2;
I_total = I_sat + DCM*parallel_axis_thm(r_cg, m_ant, I_ant)*DCM.'; % <-- need parellel axis theorem

% form R
%R = diag([I_total, I_total, I_total, I_ant, I_ant, I_ant]);
R = diag([trace(I_total), trace(I_total), trace(I_total), trace(I_ant), trace(I_ant)]);

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
DCM_LOS2ANT = simplify(DCM_LOS2SAT*DCM_SAT2ANT);

% get w_LOS2ANT in antenna frame
omega_LOS2ANT = DCM_SAT2ANT*w_sat + w_ant;


% from my diagram, 2 is roll. 1 and 3 are pitch or yaw, u choose chande
% let 1 be pitch=theta and 3=yaw=phi --> use 31 E-A
% DCM_31 = R3(e_phi) * R1(e_theta); % <-- jut for seeing relationships
e_phi_expr = acos(DCM_LOS2ANT(1,1));
e_theta_expr = acos(DCM_LOS2ANT(3,3)); % <-- this stuff feels suspect 2 me!

d_e_phi = -omega_LOS2ANT(1);  % For roll (phi)
d_e_theta = -omega_LOS2ANT(3);  % For pitch (theta)

% form constraints
phi_ant = acos(DCM_SAT2ANT(1,1));
theta_ant = acos(DCM_SAT2ANT(3,3));

max_angle = deg2rad(15); % <-- random
c = [phi_ant - max_angle; -phi_ant - max_angle; theta_ant - max_angle; -theta_ant - max_angle];

% finally form the state vector 

%% I MESSED UP AND NEED ERROR RATES STILL

dx = [
  v;
  a;
  dw_sat;
  dw_ant;
  dq_LOS2SAT;
  dq_SAT2ANT;
  d_e_phi;
  d_e_theta
];

dx = simplify(dx);


% hamiltonian
H = x.'*Q*x + u.'*R*u + lambda.'*dx + mu.'*c;
%H = simplify(H, 'IgnoreAnalyticConstraints', true, 'Steps', 50);
%disp(H)

dH_du = jacobian(H, u); 
sol = solve(dH_du.' == 0, u);
%latex(sol.u_1)

u_max = 1*ones(5,1); % <--- im assuming just |u_j| <= 1... can research
u_star = simplify(-sign(dH_du).' .* u_max);


%{

u_func = matlabFunction(u_star, 'Vars', {x1, x2, x3, x4, x5, ...
                                         lambda1, lambda2, lambda3,
                                         lambda4, lambda5});

%}