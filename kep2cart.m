%{
Convert Keplerian elements to cartesian position and velocity vectors
Inputs:
  a             - semi-major axis [m]
  e             - eccentricity
  i             - inclination [rad]
  omeg          - right ascension of ascending node [rad]
  wumbo         - argument of periapsis [rad]
  M             - mean anomaly [rad]
Outputs:
  x             - State vector in ECI frame including position and
  velocity
%}

function x_cart = kep2cart(x,m)

% Process inputs
if nargin < 2
   m = 61.6; % If no input given for m, assume zero
end

x_cart = NaN(9,length(x));

for idx = 1:length(x)
x_kep = x(:,idx);

% Retrieve orbital elements
a=x_kep(1); e=x_kep(2); i=x_kep(3);
omeg=x_kep(4); wumbo=x_kep(5); M=x_kep(6);

% Constants
G = 6.6743e-11; % Gravitational constant (m/s^)/(kg/m^2)
Me = 5.9722e24;  % Earth mass (kg)
mu = G*(Me+m);   % Gravitational constant mu

% Solve for true anomaly
E = solve_kepler(M,e); % Eccentric anomaly
cos_nu = (cos(E) - e) / (1 - e*cos(E));
sin_nu = (sqrt(1 - e^2) * sin(E)) / (1 - e*cos(E));
nu = atan2(sin_nu, cos_nu);

% Calculate semi-latus rectum
p = a*(1-e^2);

% Calculate orbital radius from focus to body
r = p/(1+e*cos(nu));

% Step 4: Position, velocity, acceleration in perifocal frame
r_perifocal = [r * cos(nu); r * sin(nu); 0];

v_perifocal = sqrt(mu / p) * [-sin(nu); e + cos(nu); 0];

a_perifocal = -mu / r^3 * r_perifocal;

% Step 5: Rotation matrices
R3_omeg  = [cos(-omeg), -sin(-omeg), 0;
            sin(-omeg),  cos(-omeg), 0;
            0,           0,          1];

R1_i  = [1, 0, 0;
         0, cos(-i), -sin(-i);
         0, sin(-i),  cos(-i)];

R3_w  = [cos(-wumbo), -sin(-wumbo), 0;
         sin(-wumbo),  cos(-wumbo), 0;
         0,                   0,    1];

Q = (R3_omeg * R1_i * R3_w)'; % total rotation

% Step 6: Rotate to ECI
r_eci = Q * r_perifocal;
v_eci = Q * v_perifocal;
a_eci = Q * a_perifocal;

x_cart(:,idx) = [r_eci; v_eci; a_eci];
end

end