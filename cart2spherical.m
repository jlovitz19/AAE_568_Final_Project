%{
Convert cartesian position and velocity vectors into spherical coordinates
Inputs:
  x_cart        - State vector in cartesian frame (6x1) [m; m/s, m/s^2]
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
a_vec = x_cart(7:9,:);

x = r_vec(1,:); y = r_vec(2,:); z = r_vec(3,:); % position
dx = v_vec(1,:); dy = v_vec(2,:); dz = v_vec(3,:); % velocity
ddx = a_vec(1); ddy = a_vec(2); ddz = a_vec(3); % acceleration

% Compute spherical coordinates
r     = vecnorm(r_vec,2);                   % radial distance
theta = atan2(y, x);                   % azimuth (longitude-like)
phi   = acos(z ./ r);                   % polar angle (from +z)

% Compute spherical velocity components
% Radial component
r_dot = dot(r_vec, v_vec) ./ r;

% Angular rates (from vector calculus in spherical coords)
theta_dot = (x.*dy - y.*dx) ./ (x.^2 + y.^2);
phi_dot   = x.*z.*dx + x.^2.*-dz + y.*(z.*dy - y.*dz);
phi_dot   = phi_dot ./ (sqrt(x.^2+y.^2) .* (x.^2+y.^2+z.^2));

% --- Second derivatives (acceleration) ---
r_ddot = (x.*ddx + y.*ddy + z.*ddz + dx.^2 + dy.^2 + dz.^2 - r_dot.^2) ./ r;

theta_ddot = (x.*ddy - y.*ddx - 2.*theta_dot.*r_dot.*(x.^2 + y.^2)) ./ (x.^2 + y.^2);

phi_ddot = ( (ddz.*r - z.*r_ddot - 2.*r_dot.*phi_dot.*r - z.*(theta_dot.^2)) ...
           - r.^2.*sin(phi).*phi_dot.^2 ) ./ (r.^2 .* cos(phi));

% Package output
x_spherical = [r; theta; phi; r_dot; theta_dot; phi_dot;...
    r_ddot; theta_ddot; phi_ddot];

end
