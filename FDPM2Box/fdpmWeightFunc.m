function [w] = fdpmWeightFunc(q, cutoff)
%fdpmWeightFunc
% Calculate weight function 
%

if exist('cutoff','var')
	q = q ./ cutoff;
end

w = (q<1) .* (1 - 6*q.^2 + 8*q.^3 - 3*q.^4);

return
end

