function [c,ceq] = AxisymEvolverConstraint(bridge, np,rp,xp)

R1 = bridge.R1;
R2 = bridge.R2;
H = bridge.H;
theta1 = bridge.theta1;
theta2 = bridge.theta2;
V = bridge.V;
X1 = bridge.X1;
X2 = bridge.X2;

r1 = rp(1);
x1 = xp(1);
r2 = rp(np);
x2 = xp(np);

c = [];
ceq = [];

%
% c(end+1:end+np) = R1 - sqrt(rp.^2 + (xp-X1).^2);
% if R2 > 0
	% c(end+1:end+np) = R2 - sqrt(rp.^2 + (xp-X2).^2);
% end
dx = (x2-x1) / (np-1);



% dr1 = rp(2)-r1;
% dl1 = sqrt(dr1^2+dx^2) * sin(theta1);
% dr2 = r2 - rp(np-1);
% dl2 = sqrt(dr2^2 + dx^2) * sin(theta2);
dl1 = 0;
dl2 = 0;

rs = rp(2:np-1);
xs = xp(2:np-1);
c(end+1:end+np-2) = R1 + dl1 - sqrt(rs.^2 + (xs-X1).^2);
if R2 > 0
	c(end+1:end+np-2) = R2 + dl2 - sqrt(rs.^2 + (xs-X2).^2);
end


%%
%%
% endpoints must on sphere
ceq(end+1) = 1 - sqrt(r1^2 + (x1-X1)^2)/R1;

% check wall
if R2 > 0
	ceq(end+1) = 1 - sqrt(r2^2 + (X2-x2)^2)/R2;
else
	ceq(end+1) = (X2 - x2)/R2;
end

% volume
% ceq(end+1) = V - AxisymEvolverVolume(bridge, np,rp,xp);
ceq(end+1) = 1 - AxisymEvolverVolume(bridge, np,rp,xp) / V;

return
end


