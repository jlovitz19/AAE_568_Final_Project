function cart_vec = spherical2cart(spherical_vec)
%Given (rho, phi, theta) find cartesian (x,y,z)
%   This vectorization should be good but if it breaks kill yourself
rho = spherical_vec(1,:);
phi = spherical_vec(2,:);
theta = spherical_vec(3,:);

x = rho.*sin(theta).*cos(phi);
y = rho.*sin(theta).*sin(phi);
z = rho.*cos(theta);

cart_vec = [x;y;z];

end

