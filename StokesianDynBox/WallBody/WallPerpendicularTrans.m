%% Resistance for sphere-wall translation in perpendicular direction
%% (Brenner, 1961)
%%

function [f] = WallPerpendicularTrans(ha,nmax)

% truncation of series
if ~exist('nmax','var')
	nmax = 500;
end

alpha = acosh(ha);
sinha = sinh(alpha);
cotha = coth(alpha);
cosecha = 1 / sinha;

ns = 1:nmax;

nc = ns .* (ns+1) ./ (2*ns-1) ./ (2*ns+3);
s1 = sinh((2*ns+1).*alpha);
s2 = (2*ns+1) .* sinh(2*alpha);
s3 = sinh((ns+0.5).*alpha);
s4 = (2*ns+1).^2 .*sinha^2;

fs = (2*s1 + s2) ./ (4*s3.^2 - s4) - 1;
fs = nc .* fs;

f = sum(fs) * sinha * 4/3;


return
end


