function [rp,xp,pres,sene] = AxisymEvolver(bridge,np,rp,xp)
% Function performing the interface optimization. 

% bridge parameters
R1 = bridge.R1;
R2 = bridge.R2;
H = bridge.H;
theta1 = bridge.theta1;
theta2 = bridge.theta2;
sigma = bridge.sigma;
X1 = bridge.X1;
X2 = bridge.X2;

% check input initial node positions
% if not provided, then generate initial guess
if isempty(rp) || isempty(xp)
	[rp,xp] = AxisymEvolverGuessInit(bridge,np);
end

% call optimizer
AxisymEvolverDriver;




return
end





