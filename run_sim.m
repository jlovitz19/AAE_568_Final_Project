clc; clear; close all;

% add attitude functions to path
run(fullfile(fileparts(mfilename('fullpath')), 'add_to_path.m'));

% Define moi of object
m_sat = 50; % kg
m_ant = 11.6; % kg
m = 61.6; % total mass - kg
l1 = .5; % m
l2 = .5; % m
I_sat = PMOI_rectangular_prism(m_sat, [l1;l1;l2]);
[I_ant,CoM_ant] = ant_inertia(0.25,0.1,0.02,m_ant);
% Use parallel axis theorem to calc total moment of inertia
I_total = I_sat + parallel_axis_thm(CoM_ant, I_ant, m_ant);
% Add parameters for bvp4c call
param.I_ant = I_ant;
param.I_total = I_total;

% Orbital constants
% Constants
G = 6.6743e-11; % Gravitational constant (m/s^)/(kg/m^2)
Me = 5.9722e24;  % Earth mass (kg)
mu = G*(Me+m);   % Gravitational constant mu
% Add parameters for bvp4c
param.mu = mu;
param.m = m;

%% Trajectory sim
% Define circular orbit IC's where orbital semi major axis is same as ISS
% Convention: 
% [semi-major axis (m); eccentricity; inclincation (rad); RAAN (rad);
% arg periapsis (rad); mean anomaly (rad)]
% Assume a = radius of earth + altitude
a = 460e3+6.37836e6; ecc = 0; inc = 0; omeg = 0; wumbo = 0; M = 0;
x0 = [a;ecc;inc;omeg;wumbo;M];

% Sim orbit returns keplerian elements
t = 0:1:100; % NOTE: SET t FOR TESTING, NORMALLY AUTOMATICALLY SET
[t,circ_orbit_kep] = sim_orbit(x0,t);
% Convert kepler --> cartesian --> spherical
% Convention: [r; theta; phi; dr; dtheta; dphi; ddr; ddtheta; ddphi]
circ_orbit_cart = kep2cart(circ_orbit_kep); % save cart coords for plotting
circ_orbit = cart2spherical(circ_orbit_cart);

%% Controller
% Set up the boundary value problem (BVP)
% Define initial guess

init_guess = NaN(1,52);  % Initial guess for state vec
init_guess(1:6) = circ_orbit(1:6,1).'; % Initial spherical coords
init_guess(7:12) = zeros(1,6); % Assume zero initial rotation

% Initial quaternions
q_LOS2SAT0 = get_initial_quat(circ_orbit(1:3,1)); % Zero rotation quat i.e. point at earth
R_31_SAT2ANT = R3(0)*R1(0);
q_SAT2ANT0 = sheppard_method(R_31_SAT2ANT); % Zero rot quat i.e. ant line with sat
init_guess(13:16) = q_LOS2SAT0; % Zero rotation quat i.e. point at earth
init_guess(17:20) = q_SAT2ANT0; % Zero rot quat i.e. ant line with sat

% Guess initial errors
% get DCM LOS2ANT 
DCM_LOS2SAT = quat2dcm(q_LOS2SAT0);
DCM_SAT2ANT = quat2dcm(q_SAT2ANT0);
DCM_LOS2ANT = DCM_LOS2SAT*DCM_SAT2ANT;
e_phi0 = acos(DCM_LOS2ANT(1,1));
e_theta0 = acos(DCM_LOS2ANT(3,3)); % <-- this stuff feels suspect 2 me!
init_guess(21:22) = [e_phi0, e_theta0]; % Initial error guess

% Guess initial extra states (x_n+1 AKA x_np1)
x_np10 = [0,0,0,0];
init_guess(23:26) = x_np10;

% Lambdas
init_guess(27:52) = zeros(1,30); % Guess initial lambdas zero

% Generate initial solution
solinit = bvpinit(t, init_guess);

% Parameters to pass to boundary value function
param.sphere0 = init_guess(1:6); % Add to parameters
param.w0 = init_guess(7:12); % angular velocities
param.q_LOS2SAT0 = q_LOS2SAT0; % quaternijawns
param.q_SAT2ANT0 = q_SAT2ANT0;
param.e_phi0 = e_phi0;       % errors
param.e_theta0 = e_theta0;
param.x_np10 = x_np10;       % extra states
param.lambda0 = init_guess(23:52); % lambdas

% Use bvp4c to solve the boundary value problem
% Wrap ODE and BC with parameters
odefun = @(x, y) bvp_ode(x, y, param);
bcfun = @(ya, yb) bvp_bcs(ya, yb, param);
opts = bvpset("AbsTol",1e-6,"RelTol", 1e-6);
sol = bvp4c(odefun, bcfun, solinit, opts);










%%
%{
Attitude component to be revised

figure
for idx = 1:4
    subplot(4,1,idx)
    plot(tv, y(:,idx))
    grid
    xlabel('time [s]', 'Interpreter','latex')
    ylabel('amplitude', 'Interpreter','latex')
    title(['$q_' num2str(idx) '$'], 'Interpreter', 'latex')
end
sgtitle('Quaternions', 'Interpreter', 'latex')

figure
for idx = 1:3
    subplot(3,1,idx)
    plot(tv, y(:,idx+4))
    grid
    xlabel('time [s]', 'Interpreter','latex')
    ylabel('km/s', 'Interpreter','latex')
    title(['$\omega_' num2str(idx) '$'], 'Interpreter', 'latex')
end    
sgtitle('Angular Velocities', 'Interpreter', 'latex')

% Plot orbit in cartesian
figure
plot3(circ_orbit_cart(1,:) / 1000 ,circ_orbit_cart(2,:)/1000,circ_orbit_cart(3,:)/1000)
title("Cartesian Plot of Simple Circular Orbit",'fontsize',12,...
    "Interpreter","Latex");
grid on;
axis equal;
xlabel("x [km]","Interpreter","Latex");
ylabel("y [km]","Interpreter","Latex");
zlabel("z [km]","Interpreter","Latex");

%}