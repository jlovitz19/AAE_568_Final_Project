%{
Convert cartesian position and velocity vectors into spherical coordinates
Inputs:
  x_cart        - State vector in cartesian frame (6x1) [m; m/s]
Outputs:
  x_spherical   - [r; theta; phi; r_dot; theta_dot; phi_dot]
                  r         - radial distance [m]
                  theta     - azimuth (longitude in xy-plane) [rad]
                  phi       - polar angle (from +z) [rad]
                  r_dot     - radial velocity [m/s]
                  theta_dot - azimuthal rate [rad/s]
                  phi_dot   - polar angle rate [rad/s]
%}
function x_spherical = cart2spherical(x_cart)

% Extract position and velocity from state vector
r_vec = x_cart(1:3,:);
v_vec = x_cart(4:6,:);

x = r_vec(1,:); y = r_vec(2,:); z = r_vec(3,:);
vx = v_vec(1,:); vy = v_vec(2,:); vz = v_vec(3,:);

% Compute spherical coordinates
r     = vecnorm(r_vec,2);                   % radial distance
theta = atan2(y, x);                   % azimuth (longitude-like)
phi   = acos(z ./ r);                   % polar angle (from +z)

% Compute spherical velocity components
% Radial component
r_dot = dot(r_vec, v_vec) ./ r;

% Angular rates (from vector calculus in spherical coords)
theta_dot = (x.*vy - y.*vx) ./ (x.^2 + y.^2);             % azimuthal rate
phi_dot   = (r.*vx.*z - x.*vz).*x + (r.*vy.*z - y.*vz).*y;   % polar rate numerator
phi_dot   = phi_dot ./ (r.^2 .* sqrt(x.^2 + y.^2));       % divide by denominator

% Package output
x_spherical = [r; theta; phi; r_dot; theta_dot; phi_dot];

end
