%% Resistance for sphere-wall rotation in perpendicular direction
%% See (Jeffery, 1916) and (Cox & Brenner, 1967)
%%

function [g] = WallPerpendicularRot(ha,nmax)

% truncation of series
if ~exist('nmax','var')
	nmax = 500;
end

alpha = acosh(ha);
sinha = sinh(alpha);
cotha = coth(alpha);
cosecha = 1 / sinha;

ns = 1:nmax;

s1 = sinh(ns.*alpha);
g = sum(s1.^(-3)) * sinha^3;

% we normalize by 6*pi*mu*a^3
g = g * 4/3;

return
end


