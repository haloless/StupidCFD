function [sene] = AxisymEvolverEnergy(bridge, np,rp,xp)

R1 = bridge.R1;
R2 = bridge.R2;
H = bridge.H;
theta1 = bridge.theta1;
theta2 = bridge.theta2;
sigma = bridge.sigma;

% contact ring radius
r1 = rp(1);
r2 = rp(np);

% embracing angle
alpha1 = asin(r1/R1);
% cap height
h1 = R1*(1-cos(alpha1));
% cap area 
scap1 = SphereCapArea(r1,h1);

% check wall case
if R2 > 0
	alpha2 = asin(r2/R2);
	h2 = R2*(1-cos(alpha2));
	scap2 = SphereCapArea(r2,h2);
else
	alpha2 = 0;
	h2 = 0;
	scap2 = pi*r2^2;
end



% rotational surface
dr = rp(2:np) - rp(1:np-1);
dx = xp(2:np) - xp(1:np-1);
% element length
dl = sqrt(dr.^2 + dx.^2);
% radius at element center
rc = 0.5 * (rp(2:np)+rp(1:np-1));
% rotational surface area
srot = sum(2*pi * rc .* dl);




% surface energy
sene = srot - cos(theta1)*scap1 - cos(theta2)*scap2;
% scale by surface tension
sene = sene * sigma;

return
end


