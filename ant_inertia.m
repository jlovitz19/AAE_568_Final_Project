%{
Calculate the antenna-body frame moment of inertia tensor of the antenna.
Optionally returns the center of mass of the antenna for use in the
parallel axis theorem.

Inputs:
    r_max: max antenna radius (m)
    d_max: max depth of the antenna from top to bottom (m)
    t_max: max thickness of antenna (at center) (m)
    m: antenna mass (kg)

Outputs
    I: 3d antenna inertia tensor in antenna body frame
    CoM: Center of mass of the antenna for use in parallel axis theorem
         ((m,m,m) coordinate)
%}

function [I,CoM] = ant_inertia(r_max,d_max,t_max,m)

a1 = (d_max - t_max)/r_max^2;
a2 = d_max/r_max^2;
b1 = t_max;

rho = (2*m)/(a1*pi*r_max^4 - a2*pi*r_max^4 + 2*b1*pi*r_max^2);

I_xx = (r_max^2*rho*pi*(a1^3*r_max^6 + 4*a1^2*b1*r_max^4 +...
    6*a1*b1^2*r_max^2 + 2*a1*r_max^4 - a2^3*r_max^6 -...
    2*a2*r_max^4 + 4*b1^3 + 3*b1*r_max^2))/12;

I_yy = (r_max^2*rho*pi*(a1^3*r_max^6 + 4*a1^2*b1*r_max^4 +...
    6*a1*b1^2*r_max^2 + 2*a1*r_max^4 - a2^3*r_max^6 -...
    2*a2*r_max^4 + 4*b1^3 + 3*b1*r_max^2))/12;

I_zz = 2*pi*((rho*(a1 - a2)*r_max^6)/6 + (b1*rho*r_max^4)/4);

I = diag([I_xx,I_yy,I_zz]);

% only give second output if requested
if nargout > 1
    z_CoM = (a1^2*r_max^4 + 3*a1*b1*r_max^2 - a2^2*r_max^4 + 3*b1^2)/...
    (6*b1 + 3*a1*r_max^2 - 3*a2*r_max^2);
    CoM = [0;0;z_CoM];
end


end