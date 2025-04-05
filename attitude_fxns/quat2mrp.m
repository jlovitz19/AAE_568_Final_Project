function [mrp] = quat2mrp(q, varargin)
%QUAT2MRP returns modified rogriguez parameters given quaternions
%   If given (q, 'shadow') it will calculate shadow set as well


mrp_curr = q(1:3) ./ (1 + q(4));

if nargin > 1 && ischar(varargin{1}) && strcmp(varargin{1}, 'shadow')
    if norm(mrp_curr) > 1
        mrp = -mrp_curr / norm(mrp_curr);  % shadow set
    else
        mrp = mrp_curr;
    end
end
mrp = mrp_curr;
end


