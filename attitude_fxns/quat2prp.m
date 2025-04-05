function [lambda,theta] = quat2prp(q)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

theta = 2*acos(q(4));
lambda = q(1:3) / sin(theta / 2);
end

