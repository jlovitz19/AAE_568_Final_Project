function [sig_dot] = mrp_kde(sig, w)
%MRP_KDE Summary of this function goes here
%   Detailed explanation goes here


sig_dot  = 1/4 * ((1 - sig.'*sig) * w + 2*cross(sig, w) + 2*w.'*sig*sig);


end

