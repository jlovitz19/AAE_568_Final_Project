function q0 = get_initial_quat(pos)
%GET_INITIAL_QUAT Summary of this function goes here
%   Detailed explanation goes here

cart_pos = spherical2cart(pos);
cart_los = cart_pos / norm(cart_pos);

default_dir = [0;0;1];

axis = cross(default_dir, cart_los);
theta = acos(dot(default_dir, cart_los));

q0 = [sin(theta/2)* axis; cos(theta/2)];


end

