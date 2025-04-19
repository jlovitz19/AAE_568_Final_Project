clc; close; clear;
G = 6.6743e-11; % Gravitational constant (m/s^)/(kg/m^2)
Me = 5.9722e24;  % Earth mass (kg)
a = 7000000; % semi major axis (m)
T = 2*pi*sqrt(a^3/(G*Me)); % orbital period (s)

[t,x] = sim_circular_orbit([a;0.5;pi/4;0;0;0],0:0.1:T,61.6);
x_cart = kep2cart(x,61.6);
x_sphere = cart2spherical(x_cart);

figure
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
