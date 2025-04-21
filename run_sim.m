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

%% Trajectory sim
% Define circular orbit IC's where orbital semi major axis is same as ISS
% Convention: 
% [semi-major axis (m); eccentricity; inclincation (rad); RAAN (rad);
% arg periapsis (rad); mean anomaly (rad)]
% Assume a = radius of earth + altitude
a = 460e3+6.37836e6; ecc = 0; inc = 0; omeg = 0; wumbo = 0; M = 0;
x0 = [a;ecc;inc;omeg;wumbo;M];

% Sim orbit returns keplerian elements
[t,circ_orbit_kep] = sim_orbit(x0);
% Convert kepler --> cartesian --> spherical
% Convention: [r; theta; phi; dr; dtheta; dphi; ddr; ddtheta; ddphi]
circ_orbit_cart = kep2cart(circ_orbit_kep); % save cart coords for plotting
circ_orbit = cart2spherical(circ_orbit_cart);

%% Controller
% Set up the boundary value problem (BVP)
% Define initial guess

init_guess = NaN(1,52);  % Initial guess for state vec
init_guess(1:6) = circ_orbit(1:6,1).'; % Initial spherical coords
init_guess(7:12) = zeros(1,6); % Assume zero initial rotation except for:
init_guess(9) = sqrt(mu)/(a^3); % w3 required for spin stabilization
init_guess(13:16) = [0,0,0,1]; % Zero rotation quat i.e. point at earth
init_guess(17:20) = [0,0,0,1]; % Zero rot quat i.e. ant line with sat
init_guess(21:22) = [0,0]; % Assume initial error and rates = 0
init_guess(23:52) = zeros(1,30);
solinit = bvpinit(t, init_guess);

% Parameters to pass to boundary value function
param.periapsis = a*(1-ecc);

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
%}

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