function [I] = PMOI_rectangular_prism(m, l)
%PMOI_RECTANGULAR_PRISM returns a 3x3 principal moment of inertia tensor
%for a rectangular prism
%   inputs:
%           m = object mass (assumed constant density)
%           l = [l1; l2; l3] lengths 
%   outputs:
%           I = inertia tensor


I = m/12 * diag([l(2)^2 + l(3)^2, l(1)^2 + l(3)^2, l(1)^2 + l(2)^2]);


end