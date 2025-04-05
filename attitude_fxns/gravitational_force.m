function Fg = gravitational_force(mu, m, R, I, order)
%GRAVITATIONAL_FORCE calculate gravitational force
%   inputs:
%       mu      :  orbital parameter
%       m       :  mass
%       R       :  radius [3x1]
%       I       :  intertia tensor
%       order   :  epsilon order of terms to include (either 1 or 2)
%   outputs:
%       Fg      : gravitational force

if order ~= [1,2]
    disp('Illegal Order')
    Fg = NaN;
else
    R_norm= norm(R);
    
    Oe2 = 3/R_norm^2 * I.'* R + 3/(2*R_norm^2) * trace(I)*R - 15/(2*R_norm^4)*(R.'*I*R)*R;
    Fg = -mu/(R_norm^3) * (m*R + Oe2*(order-1));
end

end

