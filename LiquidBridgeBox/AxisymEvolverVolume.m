function [vol] = AxisymEvolverVolume(bridge, np,rp,xp)

R1 = bridge.R1;
R2 = bridge.R2;
% H = bridge.H;
% theta1 = bridge.theta1;
% theta2 = bridge.theta2;

% contact ring
r1 = rp(1);
r2 = rp(np);

% embracing angle
alpha1 = AxisymEmbraceAngle(R1,r1);
% cap height
h1 = R1*(1-cos(alpha1));
% sphere cap volume
vcap1 = SphereCapVolume(r1,h1);

% check wall
if R2 > 0
	alpha2 = asin(r2/R2);
	h2 = R2*(1-cos(alpha2));
	vcap2 = SphereCapVolume(r2,h2);
else
	alpha2 = 0;
	h2 = 0;
	vcap2 = 0;
end

% rotational body
dx = xp(2:np) - xp(1:np-1);
% radius at element center
rc = 0.5 * (rp(2:np)+rp(1:np-1));
%
vrot = sum(pi.* rc.^2 .* dx);


% liquid volume
vol = vrot - vcap1 - vcap2;


return
end


