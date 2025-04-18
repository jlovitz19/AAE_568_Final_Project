clc; clear; close all;
syms r theta real      % cylindrical coords
syms r_max d_max t_max a1 a2 b1 % input parameters
assume([r r_max d_max t_max], 'positive')

% Paraboloid equations
z_top    = a1*r^2 + b1;
z_bottom = a2*r^2;

% Switch to cylindrical coords for integration
x = r * cos(theta);
y = r * sin(theta);
z = sym('z', 'real'); % integration variable for height

% Density (uniform)
syms rho real

% Integrate over cylindrical coords
% Limits: r ∈ [0, r_max], theta ∈ [0, 2*pi], z ∈ [z_bottom(r), z_top(r)]

% Define the integrand for total mass
mass_integrand = rho * r;
mass = int(int(int(mass_integrand, z, z_bottom, z_top), r, 0, r_max), theta, 0, 2*pi);
mass = simplify(mass);

% Define moment of inertia tensor components
% I_xx = ∫ ρ (y^2 + z^2) dV
I_xx_integrand = rho * (y^2 + z^2) * r;
I_xx = int(int(int(I_xx_integrand, z, z_bottom, z_top), r, 0, r_max), theta, 0, 2*pi);
I_xx = simplify(I_xx);

% I_yy = ∫ ρ (x^2 + z^2) dV
I_yy_integrand = rho * (x^2 + z^2) * r;
I_yy = int(int(int(I_yy_integrand, z, z_bottom, z_top), r, 0, r_max), theta, 0, 2*pi);
I_yy = simplify(I_yy);

% I_zz = ∫ ρ (x^2 + y^2) dV = ∫ ρ (r^2) dV
I_zz_integrand = rho * (r^2) * r;
I_zz = int(int(int(I_zz_integrand, z, z_bottom, z_top), r, 0, r_max), theta, 0, 2*pi);
I_zz = simplify(I_zz);

% Calculate center of mass for use in parallel axis theorem
z_integrand = rho * z * r;
z_moment = int(int(int(z_integrand, z, z_bottom, z_top), r, 0, r_max), theta, 0, 2*pi);
z_com = simplify(z_moment / mass);

clear syms r theta real r_max d_max t_max a1 a2 b1

% Calculate representative equation for density based on mas
syms m pi rho a1 a2 r_max b1

eqn = m == 2*pi*((rho*(a1 - a2)*r_max^4)/4 + (b1*rho*r_max^2)/2);

density = solve(eqn,rho);

% Output results
disp('Density in terms of mass:'); disp(density);
disp('Moment of Inertia Tensor:');
disp('I_xx ='); disp(I_xx);
disp('I_yy ='); disp(I_yy);
disp('I_zz ='); disp(I_zz);
disp('z_CoM ='); disp(z_com);
