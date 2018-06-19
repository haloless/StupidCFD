
function [k] = RecurModSphBesselK3(nmax,r,kappa)
% kappa^(n+1) * K(kappa*r)

if ~exist('kappa','var')
	kappa = 1.0;
end

% actual variable 
x = r * kappa;

f = zeros(nmax+1,2);
for n = 0:nmax
	if n == 0
		% initial value for n=0,1
		fm = 0;
		fn = 1;
		fp = x + 1;
		% derivative n=0 is -k1
		dfm = 0;
		dfn = 0;
		dfp = 1;
		% dfn = -kp;
	else
		fp = x^2 * fm + (2*n+1) * fn;
		dfp = 2*x*fm + x^2*dfm + (2*n+1)*dfn;
	end
	
	f(n+1,1) = fn;
	f(n+1,2) = dfn;
	
	fm = fn;
	fn = fp;
	dfm = dfn;
	dfn = dfp;
end

% calculate exact K function
k = zeros(nmax+1,2);
for n = 0:nmax
	fn = f(n+1,1);
	dfn = f(n+1,2);
	cn = exp(-x) / r^(n+1);
	
	% K = K(r;k)
	kn = cn * fn;
	% dK/dr
	dkn = -kappa*cn*fn - (n+1)*cn/r*fn + kappa*cn*dfn;
	
	k(n+1,1) = kn;
	k(n+1,2) = dkn;
end

if 1
	% the sqrt(pi/2) factor
	k = k .* (pi/2);
end

return
end



