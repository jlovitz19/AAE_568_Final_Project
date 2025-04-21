function y_dot = bvp_ode(t, y, param)
% ode function for boundary value problem
% 
% paramaters:
% y: vector containing state, costate, and constraint states and costates
% param: extraneous parameters used in this function, consisting of a total
% inertia matrix and an antenna inertia matrix
%
% outputs:
% y_dot: time derivative of state, costate, and constraint states and
% costates

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
    I_total = param.I_total;
    I_ant = param.I_ant;

    % cunstruct x_dot

    % IDEA FOR REPLACEMENT
    mu = param.mu;
    a_rho = -mu/rho^2;
    v = [rho_dot; rho*phi_dot; rho*th_dot*sin(phi)];

    rho_ddot = a_rho + rho*phi_dot^2 + rho*th_dot^2*sin(phi)^2;
    phi_ddot = (-2*rho_dot*phi_dot + rho*th_dot^2*sin(phi)*cos(phi))/rho;
    th_ddot = (-2*rho_dot*th_dot*sin(phi) - 2*rho*phi_dot*th_dot*cos(phi))/(rho*sin(phi));

    a = [rho_ddot; phi_ddot; th_ddot];

    % get wumbo rate of satellite body wrt inertial (note order of u_s to axis # !)
    % recall dw wants inertia tensor diagonal soo
    dw_sat = [
        (-(I_total(3,3) - I_total(2,2))*w_sat(2)*w_sat(3) + u1)/I_total(1,1);
        (-(I_total(1,1) - I_total(3,3))*w_sat(3)*w_sat(1) + u2)/I_total(2,2);
        (-(I_total(2,2) - I_total(1,1))*w_sat(1)*w_sat(2) + u3)/I_total(3,3);
    ];

    % get quanternion rate LOS2SAT
    dq_los2sat = quat_kde(q_los2sat, w_sat);

    % get antenna body w wrt sat body (note a_1 a_2 is roll. let 
    dw_ant = [
        (-(I_ant(3,3) - I_ant(2,2))*w_ant(2)*w_ant(3))/I_ant(1,1);
        (-(I_ant(1,1) - I_ant(3,3))*w_ant(3)*w_ant(1) + u4)/I_ant(2,2);
        (-(I_ant(2,2) - I_ant(1,1))*w_ant(1)*w_ant(2) + u5)/I_ant(3,3);
    ];

    % quaternion rate SAT2ANT
    dq_sat2ant = quat_kde(q_sat2ant, w_ant); 

    % get DCM LOS2ANT 
    DCM_los2sat = quat2dcm(q_los2sat);
    DCM_sat2ant = quat2dcm(q_sat2ant);
    %DCM_los2ant = simplify(DCM_los2sat*DCM_sat2ant);

    % get w_LOS2ANT
    w_los2ant = DCM_sat2ant*w_sat + w_ant;

    e_phi_dot = -w_los2ant(1);  % For roll (phi)
    e_theta_dot = -w_los2ant(3);  % For pitch (theta)

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

    lambda_dot = -get_2H2x(phi_dot, rho_dot,...
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