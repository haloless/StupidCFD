%% Generate initial guess 
function [rp,xp] = AxisymEvolverGuessInit(bridge, np)

% bridge parameters
R1 = bridge.R1;
R2 = bridge.R2;
H = bridge.H;
V = bridge.V;
theta1 = bridge.theta1;
theta2 = bridge.theta2;
sigma = bridge.sigma;
X1 = bridge.X1;
X2 = bridge.X2;

% straight cylinder solution
[rguess,aguess1,aguess2] = SolveStraightBridge(R1,R2,H,V);
% radial
rp = ones(np,1).*rguess;
% axial
if R2 > 0
	xp = linspace(X1+R1*cos(aguess1),X2-R2*cos(aguess2),np)';
else
	xp = linspace(X1+R1*cos(aguess1),X2,np)';
end


return
end

