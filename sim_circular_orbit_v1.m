function state = sim_circular_orbit_v1(h)
    % Simulates a circular orbit
    %
    % Parameters:
    % h: altitude of circular orbit above Earth's atmosphere
    %
    % Outputs:
    % state: state vector consisting of the x and y coordinates of the
    % circular orbit


    % radius of earth
    R_E = 6378;

    r = R_E + h;
    
    theta = 0:0.0001:2*pi;

    x = r*cos(theta);
    y = r*sin(theta);

    state = [x; y];
end
