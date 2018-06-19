function [c,ceq] = AxisymEvolverConstraint(bridge, np,rp,xp)
% Constraints that the optimization is subject to.

% bridge parameters
R1 = bridge.R1;
R2 = bridge.R2;
H = bridge.H;
theta1 = bridge.theta1;
theta2 = bridge.theta2;
V = bridge.V;
X1 = bridge.X1;
X2 = bridge.X2;

% two endpoints
r1 = rp(1);
x1 = xp(1);
r2 = rp(np);
x2 = xp(np);

% according to FMINCON, there are unequal (C) and equal (CEQ) constraints
c = [];
ceq = [];

%
% 0. just in case, all points should be outside the solids
%
rs = rp(2:np-1);
xs = xp(2:np-1);
c(end+1:end+np-2) = R1 - sqrt(rs.^2 + (xs-X1).^2);
if R2 > 0
	c(end+1:end+np-2) = R2 - sqrt(rs.^2 + (xs-X2).^2);
end

%
% 1. endpoints must on sphere
%
ceq(end+1) = 1 - sqrt(r1^2 + (x1-X1)^2)/R1;

% check wall
if R2 > 0
	ceq(end+1) = 1 - sqrt(r2^2 + (X2-x2)^2)/R2;
else
	ceq(end+1) = (X2 - x2)/R2;
end

%
% 2. volume must be conserved
%
ceq(end+1) = 1 - AxisymEvolverVolume(bridge, np,rp,xp) / V;

return
end


