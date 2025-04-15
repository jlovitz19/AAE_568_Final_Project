%{
Solves Kepler's Equation for Eccentric Anomaly
  Inputs:
    M         - Mean anomaly (in radians)
    e         - Orbital eccentricity (0 < e < 1)
    tol       - Tolerance for convergence (default: 1e-10)
    max_iter  - Maximum iterations (default: 100)
  
  Output:
    E         - Eccentric anomaly (in radians)
%}

function E = solve_kepler(M,e,tol,max_iter)

% Set default tolerance and max iterations if inputs not given
if nargin < 3
    tol = 1e-3;
end
if nargin < 4
    max_iter = 10;
end

% Normalize M to the range [0, 2*pi)
M = mod(M, 2*pi);

% Initial guess for E
if e < 0.8
    E = M;
else
    E = pi;
end

% Newton-Raphson iteration
for i = 1:max_iter
    f = E - e*sin(E) - M;
    f_prime = 1 - e*cos(E);
    delta = -f / f_prime;
    E = E + delta;

    if abs(delta) < tol
        return;
    end
end

warning('solve_kepler: Did not converge within %d iterations', max_iter);
end
