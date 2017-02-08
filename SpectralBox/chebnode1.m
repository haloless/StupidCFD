
function [ xs ] = chebnode1(n,xa,xb)
% Chebyshev points of the 1st kind 
% xk = cos(pi*k/n), 0<=k<=n
% It contains two boundary nodes at endpoints +1 and -1.
%
% NOTE this is different from Chebyshev points of the 2nd kind
%
if ~exist('xa','var')
	xa = -1.0;
end
if ~exist('xb','var')
	xb = 1.0;
end

if n == 0
	xs = 0;
else
	ks = 0:n-1;
	% definition by cosine
	% xs = cos(pi/(n-1) .* ks);
	% sine version, numerically more symmetric!
	xs = sin(pi/(2*n) * (n-1-ks*2));
	% xs = (xb-xa)/2 .* (xs+1) + xa;
	xs = xs';
end



return
end

