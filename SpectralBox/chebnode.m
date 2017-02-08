
function [ xs ] = chebnode(n,xa,xb)
% Chebyshev points of the 2nd kind
% (or Chebyshev extreme point)
% (or Gauss-Chebyshev point)
%
% xk = cos(pi*k/n), 0<=k<=n
% It contains two boundary nodes at endpoints +1 and -1.
%
% NOTE this is different from Chebyshev points of the 1st kind
%
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

