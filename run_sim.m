clc; clear; close all;

% add attitude functions to path
run(fullfile(fileparts(mfilename('fullpath')), 'add_to_path.m'));


% define moi of object
m = 250; %kg
l1 = .5; %m
l2 = .75;
I = PMOI_rectangular_prism(m, [l1;l1;l2] / 1000);

% Define circular orbit IC's where orbital semi major axis is same as ISS
% Convention: 
%     [semi-major axis (m); eccentricity; inclincation (rad); RAAN (rad);
%      arg periapsis (rad); mean anomaly (rad)]
a = 460e3+6.37836e6; % semi major axis --> radius of Earth plus alt.
x0 = [a;0;0;0;0;0];

% Sim orbit returns keplerian elements
[t,circ_orbit_kep] = sim_orbit(x0);
% Convert kepler --> cartesian --> spherical
% Convention: [r; theta; phi; dr; dtheta; dphi; ddr; ddtheta; ddphi]
circ_orbit_cart = kep2cart(circ_orbit_kep); % save cart coords for plotting
circ_orbit = cart2spherical(circ_orbit_cart);

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