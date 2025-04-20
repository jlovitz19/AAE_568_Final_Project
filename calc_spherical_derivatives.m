clc; clear; close;
% Calculate the symbolic representations of the spherical coordinates and
% their derivatives, expressed in cartesian coordinates for coordinate
% conversion

syms x(t) y(t) z(t) t r(t) th(t) phi(t) dx dy dz ddx ddy ddz

r(t) = sqrt(x(t)^2+y(t)^2+z(t)^2);
th(t) = atan(y(t)/x(t));
phi(t) = acos(z/(sqrt(x(t)^2+y(t)^2+z(t)^2)));

dr = diff(r(t),t);
ddr = diff(diff(r(t),t),t);

dphi = diff(phi(t),t);
ddphi = diff(diff(phi(t),t),t);

dth = diff(th(t),t);
ddth = diff(diff(th(t),t),t);
