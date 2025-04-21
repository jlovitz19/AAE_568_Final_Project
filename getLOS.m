function [LOS] = getLOS(pos_spherical)
% [r, phi, theta]
% Gets line of sight toward Earth vector for circular coordinates
% Ensure vectors are input as num states x N

LOS = -pos_spherical(1:3,:) ./ vecnorm(pos_spherical(1:3,:));
end

