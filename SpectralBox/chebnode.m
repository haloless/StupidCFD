
function [ xs ] = chebnode(n,xa,xb)

if ~exist('xa','var')
	xa = -1.0;
end
if ~exist('xb','var')
	xb = 1.0;
end

if n == 0
	xs = 1.0;
else
	xs = cos(pi/n .* (0:n))';
	xs = (xb-xa)/2 .* (xs+1) + xa;
end



return
end

