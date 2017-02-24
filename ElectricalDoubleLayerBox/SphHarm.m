
function [Ynm] = SphHarm(n,m,theta,phi)

szt = size(theta);
szp = size(phi);

theta = theta(:);
phi = phi(:);

absm = abs(m);
mu = cos(theta);

lp = legendre(n,mu);
lp = lp(absm+1,:)';

N = sqrt(factorial(n-absm)/factorial(n+absm));

em = exp(1i * m * phi);

Ynm = N .* lp .* em;

if 1
	Ynm = Ynm * (-1)^m;
end

if 1
	Ynm = Ynm * sqrt((2*n+1)/(4*pi));
end

if numel(theta) > 1
	Ynm = reshape(Ynm,szt);
elseif numel(phi) > 1
	Ynm = reshape(Ynm,szp);
end

return
end

