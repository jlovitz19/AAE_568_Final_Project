clc; clear; close all;

script_path = matlab.desktop.editor.getActiveFilename;
run(fullfile(fileparts(mfilename('fullpath')), 'add_to_path.m'));

syms rho theta phi drho dtheta dphi ddrho ddtheta ddphi omega_sat_1 omega_sat_2 omega_sat_3 omega_ant_1 omega_ant_2 omega_ant_3 
syms q_sat_1 q_sat_2 q_sat_3 q_sat_4 q_ant_1 q_ant_2 q_ant_3 q_ant_4 e_theta e_phi real
syms lambda_1 lambda_2 lambda_3 lambda_4 lambda_5 lambda_6 lambda_7 lambda_8 lambda_9 lambda_10 lambda_11 lambda_12 lambda_13 lambda_14 lambda_15 lambda_16 lambda_17 lambda_18 lambda_19 lambda_20 lambda_21 lambda_22 lambda_23 lambda_24 lambda_25 lambda_26
syms mu_1 mu_2 mu_3 mu_4 mu_5 mu_6 mu_7 mu_8 mu_9 mu_10 u_1 u_2 u_3 u_4 u_5 
syms x_np1_1 x_np1_2 x_np1_3 x_np1_4
assume(q_sat_1^2 + q_sat_2^2 + q_sat_3^2 + q_sat_4^2 == 1)
assume(q_ant_1^2 + q_ant_2^2 + q_ant_3^2 + q_ant_4^2 == 1)

l = .5;
m_sat = 50;
r_max = 0.25;
d_max = 0.05;
t_max = 0.01;
m_ant = 11.6;
u_ant_max = 0.5; % Nm... arbitrary.. roughly 0.1-2 Nm range
u_sat_max = 0.5; % Nm... arbitrary (and different from ant_max) roughly 0.05-0.5 Nm
error_weight = 10; % arbitrary.. to be tuned

% form state
w_sat = [omega_sat_1, omega_sat_2, omega_sat_3].';
w_ant = [omega_ant_1, omega_ant_2, omega_ant_3].';
q_LOS2SAT = [q_sat_1, q_sat_2, q_sat_3, q_sat_4].';
q_SAT2ANT = [q_ant_1, q_ant_2, q_ant_3, q_ant_4].';
err = [e_phi, e_theta].';

%% ORDER IS RHO, PHI, THETA %%

% erm are dphi dro and dtheta supposed ?
x_np1 = [x_np1_1; x_np1_2; x_np1_3; x_np1_4];
x = [rho; phi; theta; drho; dphi; dtheta; w_sat; w_ant; q_LOS2SAT; q_SAT2ANT; err; x_np1]; 

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

%{
I_tot = matlabFunction([I_total(1, 1); I_total(2, 2); I_total(3, 3)], 'File', 'get_totalInertia.m');
I_a = matlabFunction([I_ant(1, 1); I_ant(2, 2); I_ant(3, 3)], 'File', 'get_antennaInertia.m');
%}
% form R
%R = diag([I_total, I_total, I_total, I_ant, I_ant, I_ant]);
%R = diag([trace(I_total), trace(I_total), trace(I_total), trace(I_ant), trace(I_ant)]);
R_ = diag([I_total(1,1), I_total(2,2), I_total(3,3)]);
R = blkdiag(R_, I_ant(1:2, 1:2)); % davis correct me if im wrong but I(3,3) would be roll?

% form costate
lambda = [lambda_1 lambda_2 lambda_3 lambda_4 lambda_5 lambda_6 lambda_7 lambda_8 lambda_9 lambda_10 lambda_11 lambda_12 lambda_13 lambda_14 lambda_15 lambda_16 lambda_17 lambda_18 lambda_19 lambda_20 lambda_21 lambda_22 lambda_23 lambda_24 lambda_25 lambda_26].';
mu = [mu_1 mu_2 mu_3 mu_4 mu_5 mu_6 mu_7 mu_8 mu_9 mu_10].';

%% dynamics !

G = 6.6743e-11; % Gravitational constant (m/s^)/(kg/m^2)
Me = 5.9722e24;  % Earth mass (kg)
GM = G*(Me+61.6);   % Gravitational constant mu

% get vel and acc
%v = [drho; rho*dphi; rho*dtheta*sin(phi)];
v = [drho; dphi; dtheta];

a_rho = -GM/rho^2;
ddrho = a_rho + rho*dphi^2 + rho*dtheta^2*sin(phi)^2;
ddphi = (-2*drho*dphi + rho*dtheta^2*sin(phi)*cos(phi))/rho;
ddtheta = (-2*drho*dtheta*sin(phi) - 2*rho*dphi*dtheta*cos(phi))/(rho*sin(phi));

a = [ddrho; ddphi; ddtheta];

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
c_ant_angle = [phi_ant - max_angle; -phi_ant - max_angle; theta_ant - max_angle; -theta_ant - max_angle];

c_max_torque = [
    u_1 - u_sat_max; -u_1 - u_sat_max;
    u_2 - u_sat_max; -u_2 - u_sat_max;
    u_3 - u_sat_max; -u_3 - u_sat_max;
    u_4 - u_ant_max; -u_4 - u_ant_max
    u_5 - u_ant_max; -u_5 - u_ant_max
    ];
c = c_max_torque; %[c_ant_angle; c_max_torque];
% finally form the state vector 
dx = [
  v;
  a;
  dw_sat;
  dw_ant;
  dq_LOS2SAT;
  dq_SAT2ANT;
  d_e_phi;
  d_e_theta;
  c_ant_angle
];

dx = simplify(dx);


% hamiltonian
H = x.'*Q*x + u.'*R*u + lambda.'*dx + mu.'*c;
%H = simplify(H, 'IgnoreAnalyticConstraints', true, 'Steps', 50);
%disp(H)

dH_du = jacobian(H, u);
%dH_du = simplify(dH_du, 'IgnoreAnalyticConstraints', true, 'Steps', 50);
sol = solve(dH_du.' == 0, u);
%sol = struct2array(sol);
%sol = simplify(sol, 'IgnoreAnalyticConstraints', true, 'Steps', 50);

mu_soln = solve(dH_du.' == 0, mu);
names = fieldnames(mu_soln);
for idx = 1:numel(names)
    matlabFunction(mu_soln.(names{idx}), 'File', ['get_', names{idx}]);
end


u_1_func = matlabFunction(sol.u_1, 'File', 'get_u1.m');
u_2_func = matlabFunction(sol.u_2, 'File', 'get_u2.m');
u_3_func = matlabFunction(sol.u_3, 'File', 'get_u3.m');
u_4_func = matlabFunction(sol.u_4, 'File', 'get_u4.m');
u_5_func = matlabFunction(sol.u_5, 'File', 'get_u5.m');

%latex(sol.u_1)
%{
u_max = 1*ones(5,1); % <--- im assuming just |u_j| <= 1... can research
u_star = simplify(-sign(dH_du).' .* u_max);
%}
%{
u_bangbang = matlabFunction(u_star, 'Vars', {rho theta phi ...
drho dtheta dphi ddrho ddtheta ddphi omega_sat_1 omega_sat_2 omega_sat_3 ...
omega_ant_1 omega_ant_2 omega_ant_3 q_sat_1 q_sat_2 q_sat_3 q_sat_4  ... 
q_ant_1 q_ant_2 q_ant_3 q_ant_4 e_theta e_phi error_weight
});
%}
%u_bangbang = matlabFunction(u_star);


%{

u_func = matlabFunction(u_star, 'Vars', {x1, x2, x3, x4, x5, ...
                                         lambda1, lambda2, lambda3,
                                         lambda4, lambda5});

%}

% lambda dot = 2H/2x
del_H_del_x = jacobian(H, x);
get_del_h_del_x = matlabFunction(del_H_del_x, 'File', 'get_2H2x');
lambda_dot = -del_H_del_x;

% 2H/2u = 0
del_H_del_u = jacobian(H, u) == zeros(1, 5);
u1 = solve(del_H_del_u(1, 1), u_1);
u2 = solve(del_H_del_u(1, 2), u_2);
u3 = solve(del_H_del_u(1, 3), u_3);
u4 = solve(del_H_del_u(1, 4), u_4);
u5 = solve(del_H_del_u(1, 5), u_5);


%%% end of E-L equations %%%

%%% start of 2-pt bvp %%%



%{

u_func = matlabFunction(u_star, 'Vars', {x1, x2, x3, x4, x5, ...
                                         lambda1, lambda2, lambda3,
                                         lambda4, lambda5});

%}
function y_dot = bvp_ode(y, p)
    % antenna limits (can be changed)
    gamma_1 = 15*pi/180;
    gamma_2 = 15*pi/180;
    gamma_3 = 0.5;
    gamma_4 = 0.5;

    % states
    rho = y(1);
    phi = y(2);
    theta = y(3);

    rho_dot = y(4);
    phi_dot = y(5);
    th_dot = y(6);

    w_sat = y(7:9);
    w_ant = y(10:12);

    q_los2sat = y(13:16);
    q_sat2ant = y(17:20);

    q_sat_1 = q_los2sat(1);
    q_sat_2 = q_los2sat(2);
    q_sat_3 = q_los2sat(3);
    q_sat_4 = q_los2sat(4);

    q_ant_1 = q_sat2ant(1);
    q_ant_2 = q_sat2ant(2);
    q_ant_3 = q_sat2ant(3);
    q_ant_4 = q_sat2ant(4);

    e_phi = y(21);
    e_theta = y(22);

    % costates
    lambda_1 = y(27);
    lambda_2 = y(28);
    lambda_3 = y(29);
    lambda_4 = y(30);
    lambda_5 = y(31);
    lambda_6 = y(32);
    lambda_7 = y(33);
    lambda_8 = y(34);
    lambda_9 = y(35);
    lambda_10 = y(36);
    lambda_11 = y(37);
    lambda_12 = y(38);
    lambda_13 = y(39);
    lambda_14 = y(40);
    lambda_15 = y(41);
    lambda_16 = y(42);
    lambda_17 = y(43);
    lambda_18 = y(44);
    lambda_19 = y(45);
    lambda_20 = y(46);
    lambda_21 = y(47);
    lambda_22 = y(48);

    % state constraint states and costates
    x_np1 = y(23:26);
    lambda_23 = y(49);
    lambda_24 = y(50);
    lambda_25 = y(51);
    lambda_26 = y(52);

    % unconstrained inputs
    u1 = get_u1(lambda_7, 0, 0, q_ant_1, q_ant_2, q_ant_3, q_ant_4);
    u2 = get_u2(lambda_8, 0, 0, q_ant_1, q_ant_2, q_ant_3, q_ant_4);
    u3 = get_u3(lambda_9, 0, 0, q_ant_1, q_ant_2, q_ant_3, q_ant_4);
    u4 = get_u4(lambda_11, 0, 0);
    u5 = get_u5(lambda_12, 0, 0);

    % get inertia
    I_total = p.I_total;
    I_ant = p.I_ant;

    % cunstruct x_dot

    % get vel and acc
    v = [rho_dot; rho*phi_dot; rho*th_d*sin(phi)];
    a = [
        (rho_ddot - rho*phi_d^2 - rho*th_d^2*sin(phi)^2);
        rho*phi_ddot + 2*rho_dot*phi_dot - rho*th_dot^2*sin(phi)*cos(phi);
        rho*th_ddot*sin(phi) + 2*rho_dot*th_dot*sin(phi) + 2*rho*phi_dot*th_dot*cos(phi)
    ];

    % get wumbo rate of satellite body wrt inertial (note order of u_s to axis # !)
    % recall dw wants inertia tensor diagonal soo
    dw_sat = [
        (-(I_total(3,3) - I_total(2,2))*w_sat(2)*w_sat(3) + u1)/I_total(1,1);
        (-(I_total(1,1) - I_total(3,3))*w_sat(3)*w_sat(1) + u_2)/I_total(2,2);
        (-(I_total(2,2) - I_total(1,1))*w_sat(1)*w_sat(2) + u_3)/I_total(3,3);
    ];

    % get quanternion rate LOS2SAT
    dq_los2sat = quat_kde(q_los2sat, w_sat);

    % get antenna body w wrt sat body (note a_1 a_2 is roll. let 
    dw_ant = [
        (-(I_ant(3,3) - I_ant(2,2))*w_ant(2)*w_ant(3))/I_ant(1,1);
        (-(I_ant(1,1) - I_ant(3,3))*w_ant(3)*w_ant(1) + u_4)/I_ant(2,2);
        (-(I_ant(2,2) - I_ant(1,1))*w_ant(1)*w_ant(2) + u_5)/I_ant(3,3);
    ];

    % quaternion rate SAT2ANT
    dq_sat2ant = quat_kde(q_sat2ant, w_ant); 

    % get DCM LOS2ANT 
    DCM_los2sat = quat2dcm(q_los2sat);
    DCM_sat2ant = quat2dcm(q_sat2and);
    %DCM_los2ant = simplify(DCM_los2sat*DCM_sat2ant);

    % get w_LOS2ANT
    w_los2ant = DCM_sat2ant*w_sat + w_ant;

    e_phi_dot = -w_los2ant(1);  % For roll (phi)
    e_theta_dot = -w_los2and(3);  % For pitch (theta)

    % form constraints
    phi_ant = acos(DCM_sat2ant(1,1));
    theta_ant = acos(DCM_sat2ant(3,3));

    C = [u1-gamma_3; -u1-gamma_3; u2-gamma_3; -u2-gamma_3;...
        u3-gamma_3; -u3-gamma_3; u4-gamma_4; -u4-gamma_4;...
        u5-gamma_4; -u5-gamma_4];

    %{
    u1 = get_u1(lambda_7, mu_1, mu_2, q_ant_1, q_ant_2, q_ant_3, q_ant_4);
    u2 = get_u2(lambda_8, mu_3, mu_4, q_ant_1, q_ant_2, q_ant_3, q_ant_4);
    u3 = get_u3(lambda_9, mu_5, mu_6, q_ant_1, q_ant_2, q_ant_3, q_ant_4);
    u4 = get_u4(lambda_11, mu_7, mu_8);
    u5 = get_u5(lambda_12, mu_9, mu_10);
    %}

    mu = zeros(10, 1);

    if C(1) >= 0
        u1 = gamma_3;
        mu(1) = get_mu_1(lambda_7, q_ant_1,...
            q_ant_2, q_ant_3, q_ant_4, u1);
    elseif C(2) >= 0
        u1 = -gamma_3;
        mu(2) = get_mu_2;
    elseif C(3) >= 0
        u2 = gamma_3;
        mu(3) = get_mu_3(lambda_8, q_ant_1,...
            q_ant_2, q_ant_3, q_ant_4, u2);
    elseif C(4) >= 0
        u2 = -gamma_3;
        mu(4) = get_mu_4;
    elseif C(5) >= 0
        u3 = gamma_3;
        mu(5) = get_mu_5(lambda_9, q_ant_1,...
            q_ant_2, q_ant_3, q_ant_4, u3);
    elseif C(6) >= 0
        u3 = -gamma_3;
        mu(6) = get_mu_6;
    elseif C(7) >= 0
        u4 = gamma_4;
        mu(7) = get_mu_7(lambda_11, u4);
    elseif C(8) >= 0
        u4 = -gamma_4;
        mu(8) = get_mu_8;
    elseif C(9) >= 0
        u5 = gamma_4;
        mu(9) = get_mu_9(lambda_12, u5);
    elseif C(10) >= 0
        u5 = -gamma_4;
        mu(10) = get_mu_10;
    end
 
    % heaveside step functions for state constraint
    if phi_ant - gamma_1 >= 0
        I_g1 = 0;
    else
        I_g1 = 1;
    end

    if -phi_ant - gamma_1 >= 0
        I_g2 = 0;
    else
        I_g2 = 1;
    end

    if theta_ant - gamma_2 >= 0
        I_g3 = 0;
    else
        I_g3 = 1;
    end

    if -theta_ant - gamma_2 >= 0
        I_g4 = 0;
    else
        I_g4 = 1;
    end

    % state constrained x_dot
    x_dot_np1 = [(phi_ant - gamma_1)^2*I_g1;...
        (-phi_ant - gamma_1)^2*I_g2;...
        (theta_ant - gamma_2)^2*I_g3;...
        (-theta_ant - gamma_2)^2*I_g4];

    % construct x_dot
    x_dot = [v; a; dw_sat; dw_ant; dq_los2sat;...
        dq_sat2ant; e_phi_dot; e_theta_dot; x_dot_np1];

    lambda_dot = -get_2H2x(phi_ddot, th_ddot, phi_dot, rho_dot,...
        th_dot, lambda_1, lambda_2, lambda_3, lambda_4, lambda_5,...
        lambda_6, lambda_7, lambda_8, lambda_9, lambda_10,...
        lambda_11, lambda_13, lambda_14, lambda_15, lambda_16,...
        lambda_17, lambda_18, lambda_19, lambda_20, lambda_21,...
        lambda_22, lambda_23, lambda_24, lambda_25, lambda_26,...
        w_ant(1), w_ant(2), w_ant(3), w_sat(1), w_sat(2), w_sat(3),...
        phi, q_ant_1, q_ant_2, q_ant_3, q_ant_4, q_sat_1, q_sat_2,...
        q_sat_3, q_sat_4, rho, u1, u2, u3, x_np1(3), x_np1(4));

    y_dot = [x_dot; lambda_dot]';

end

function bcs = bvp_bcs(yi, yf)
    R = 6000000; % can be changed based on rho IC

    bcs = [yi(1)-R, yi(2), yi(3), yi(4), yi(5), yi(6),...
        yi(7), yi(8), yi(9), yi(10), yi(11), yi(12), yi(13),...
        yi(14), yi(15), yi(16), yi(17), yi(18), yi(19), yi(20),...
        yi(21), yi(22), yf(1)-yi(1), yf(2)-yi(2), yf(3)-yi(3),...
        yf(4)-yi(4), yf(5)-yi(5), yf(6)-yi(6),...
        yf(7), yf(8), yf(9), yf(10), yf(11), yf(12), yf(13),...
        yf(14), yf(15), yf(16), yf(17), yf(18), yf(19), yf(20), yf(21),...
        yf(22), yf(23), yi(23), yf(24), yi(24), yf(25), yi(25),...
        yf(26), yi(26)];
end