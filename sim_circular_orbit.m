% Sim circular orbit
%{
Inputs
Initial state "x0" representing classical orbital elements:
    a0: initial semi-major axis (km)
    e0: initial eccentricty
    i0: initial inclination (rad)
    omeg0: initial right ascension (rad)
    wumbo0: initial argument of periapsis (rad)
    M0: initila mean anomaly (rad)
tspan: time span over which to integrate
m: mass of the satellite
a_d: control input (vector)
%}

function [t,x] = sim_circular_orbit(x0,tspan,m,a_d)

% Process inputs
if nargin < 3
   m = 0; % If no input given for m, assume zero
end
if nargin < 4
    a_d = zeros(3,1); % If no control input given --> assume no control
end

options = odeset(RelTol=1e-8, AbsTol=1e-8);
[t,x] = ode45(@(t,x) Pert_2B_EOM(x,m,a_d),tspan,x0,options);
t = t.';
x = x.';

end