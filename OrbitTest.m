clc; close; clear;
G = 6.6743e-11; % Gravitational constant (m/s^)/(kg/m^2)
Me = 5.9722e24;  % Earth mass (kg)
% Radius of the Earth (can be 1 for normalized sphere)
R = 6.37836e6;
a = 460e3+R; % semi major axis (m)
T = 2*pi*sqrt(a^3/(G*Me)); % orbital period (s)

[t,x] = sim_circular_orbit([a;0.00054;deg2rad(51.6360);deg2rad(15.8423);...
    deg2rad(55.3016);deg2rad(304.8477)],0:0.1:T,61.6);
x_cart = kep2cart(x,61.6);
x_sphere = cart2spherical(x_cart);

% Create Earth for plotting with orbit
% Load an Earth texture image
earth = imread('2k_earth_daymap.jpg');

% Create spherical coordinates
[lon,lat] = meshgrid(linspace(-pi, pi, size(earth,2)),...
    linspace(-pi/2, pi/2, size(earth,1)));

% Convert to Cartesian coordinates
xe = R * cos(lat) .* cos(lon);
ye = R * cos(lat) .* sin(lon);
ze = R * sin(lat);

% Plot orbit and Earth
figure
surf(xe,ye,ze,'EdgeColor','none')
hold on;
colormap('gray');  % Make sure the surface uses the correct color map
shading interp;    % Smooth color transitions
surface(xe,ye,ze,'FaceColor', 'texturemap', 'CData', flipud(earth),...
    'EdgeColor', 'none');
plot3(x_cart(1,:),x_cart(2,:),x_cart(3,:))

axis equal
xlabel("x")
ylabel("y")
zlabel("z")
grid on;

figure
names = ["$\rho$","$\theta$","$\phi$"];
for idx = 1:3
subplot(3,1,idx)
plot(t,x_sphere(idx,:));
xlabel("time (s)",'fontsize',20,"Interpreter","Latex");
ylabel(names(idx),'fontsize',20,"Interpreter","Latex")
grid on;
end
