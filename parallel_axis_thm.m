function I = parallel_axis_thm(r_CG, I_C, m)
    % calculates I of a body using parallel axis theorem
    %
    % Parameters:
    % r_CG: distance of body from center of mass of system
    % I_C: mass of individual body
    % m: mass of body
    %
    % Outputs:
    % T: new inertia matrix post-PAT

    if nargin < 3
        m = 11.6;
    end

    x_G = r_CG(1);
    y_G = r_CG(2);
    z_G = r_CG(3);

    mat_comp = [y_G^2+z_G^2 -x_G*y_G -x_G*z_G;...
        -x_G*y_G x_G^2+z_G^2 -y_G*z_G;...
        -x_G*z_G -y_G*z_G x_G^2+y_G^2];

    I = I_C + m*mat_comp;
end
