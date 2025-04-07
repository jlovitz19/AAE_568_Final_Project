function x_dot = sim_KDE_DDE(x, J)
    % Simulates the kinematic differential equation and the dynamic 
    % differential equation simultaneously and returns the attitude
    % represented using quaternions along with the object's angular
    % velocity
    %
    % Parameters:
    % R: orbital radius from the center of Earth
    %
    % Outputs:
    % q_dot: quaternion attitude representation time history
    % wumbo_dot: angular velocity vector of spacecraft (rad/s)

    % retrieve states
    q1 = x(1);
    q2 = x(2);
    q3 = x(3);
    q4 = x(4);

    w1 = x(5);
    w2 = x(6);
    w3 = x(7);

    % retrieve moments of inertia
    J1 = J(1);
    J2 = J(2);
    J3 = J(3);

    % KDE
    q_dot = .5*[q4 -q3 q2 q1; q3 q4 -q1 q2; -q2 q1 q4 q3;...
        -q1 -q2 -q3 q4]*[w1; w2; w3; 0];

    % DDE
    wumbo_dot = [(J2-J3)*w2*w3/J1; (J3-J1)*w1*w3/J2; (J1-J2)*w1*w2/J3];

    % combine states
    x_dot = [q_dot; wumbo_dot];
end

