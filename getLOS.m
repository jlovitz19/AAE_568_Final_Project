function [LOS] = getLOS(pos_spherical)
% [r, theta, phi]
%   Detailed explanation goes here
LOS = -pos_spherical(1,:) / norm(pos_spherical(1,:));
end

