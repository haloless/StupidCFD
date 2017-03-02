function [Ynm] = SphHarmY(n,m,theta,phi)
% Spherical Harmonics Y_n^m(theta,phi)
%

szt = size(theta);
szp = size(phi);

theta = theta(:);
phi = phi(:);

absm = abs(m);
mu = cos(theta);

% associated Legendre function
pnm = legendre(n,mu);
pnm = pnm(absm+1,:)';

N = sqrt(factorial(n-absm)/factorial(n+absm));

em = exp(1i * m * phi);

Ynm = N .* pnm .* em;

if 1
	% the phase factor
	Ynm = Ynm * (-1)^m;
end

if 0
	% normalization factor
	Ynm = Ynm * sqrt((2*n+1)/(4*pi));
end

if numel(theta) > 1
	Ynm = reshape(Ynm,szt);
elseif numel(phi) > 1
	Ynm = reshape(Ynm,szp);
end

return
end

