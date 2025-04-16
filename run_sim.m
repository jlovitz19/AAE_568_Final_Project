clc; clear; close all;

% add attitude functions to path
run(fullfile(fileparts(mfilename('fullpath')), 'add_to_path.m'));

% define physical parameters
h = 1000; %km

% define moi of object
m = 250; %kg
l1 = .5; %m
l2 = .75;
I = PMOI_rectangular_prism(m, [l1;l1;l2] / 1000);

% define time
tv = linspace(0, 2000, 1000);

% define IC (arbitrarily)
% x is quaternions, wumbos
x0 = [0;0;0;1; .01; -.01; .1];


% solve ode
opts = odeset('RelTol',1e-6, 'AbsTol',1e-6);
[~, y] = ode45(@(t,y) sim_KDE_DDE(y, diag(I)), tv, x0, opts);

circ_orbit = sim_circular_orbit_v1(h);


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

figure
plot(circ_orbit(1,:), circ_orbit(2,:))
xlabel('km', 'Interpreter', 'latex')
ylabel('km', 'Interpreter', 'latex')
grid
title('Circular Orbit', 'Interpreter', 'latex')
axis equal