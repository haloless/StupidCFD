function [bridge] = MakeBridge(R1,R2,H,theta1,theta2,V,sigma)
% Create a bridge struct holding given parameters.

if ~exist('sigma','var')
	sigma = 1.0;
end

%
bridge = struct();

bridge.R1 = R1;
bridge.R2 = R2;
bridge.H = H;
bridge.theta1 = theta1;
bridge.theta2 = theta2;
bridge.V = V;
bridge.sigma = sigma;

% two sphere centers
bridge.X1 = 0;
if R2 > 0
	bridge.X2 = R1 + H + R2;
else
	bridge.X2 = R1 + H;
end


return
end


