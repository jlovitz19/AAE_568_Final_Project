%{
Perturbed two body dynamics via Gauss EoM's on classical orbital elements

Inputs
State "x" representing classical orbital elements:
    a: semi-major axis (m)
    e: eccentricty
    i: inclination (rad)
    omeg: right ascension (rad)
    wumbo: argument of periapsis (rad)
    M: mean anomaly (rad)
m: mass of the satellite
a_d: control input

Outputs
x_dot: the time derivative of the state

Assumptions:
 - Earth centered orbit
 - Using radians

Note: must add "solve_kepler.m" to path in main runner script
%}

function x_dot = Pert_2B_EOM(x,m,a_d)

% Retrieve orbital elements
a=x(1); e=x(2); i=x(3); omeg=x(4); wumbo=x(5); M=x(6);

% Process inputs
if nargin < 2
   m = 61.6; % If no input given for m, assume 61.6
end
if nargin < 3
   a = zeros(3,1); % If no input given for m, assume zero
end
    
% Constants
G = 6.6743e-11; % Gravitational constant (km/s^)/(kg/m^2)
Me = 5.9722e24;  % Earth mass (kg)

% Mean motion
n = sqrt(G*(Me+m)/a^3);
f0 = [0,0,0,0,0,n]';

x_dot = f0;

% Check to see if control input is even applicable
if norm(a_d) ~= 0

% Calculate angular momentum
mu = Me*m/(Me+m); % reduced mass
h = mu*sqrt(G*Me*a*(1-e^2));

% Solve for true anomaly
E = solve_kepler(M,e); % Eccentric anomaly
nu = 2*atan2(sqrt((1+e)*tan(E/2)),sqrt(1-e));

% Semi-minor axis
b = a*sqrt(1-e^2);

% Calculate semi-latus rectum
p = a*(1-e^2);

% Calculate orbital radius from focus to body
r = p/(1+e*cos(nu));

% Define B matrix
B = 1/h * [2*a^2*e*sin(nu),  2*a^2*p/r,  0;
           p*sin(nu),  (p+r)*cos(nu)+r*e,  0;
           0,  0,  r*cos(nu+wumbo);
           0,  0,  r*sin(nu+wumbo/sin(i));
           -p*cos(nu)/e,  (p+r)*sin(nu)/e,  -r*sin(nu+wumbo)/tan(i);
           b*p*cos(nu)/(a*e) - 2*b*r/a,  -b*(p+r)*sin(nu)/(a*e),  0];

% Calculate time derivative of state
x_dot = x_dot + B*a_d;

end

end