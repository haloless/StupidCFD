
clear all;

sigma = 1.0;
theta = deg2rad(15);
theta1 = theta;
theta2 = theta;
R1 = 1.0;
R2 = 1.0;
% H = 0.2;
V = SphereVolume(R1)*0.005;

dist = linspace(0,0.2,21);

if 1
	for H = dist
		
	end
end

bridge = MakeBridge(R1,R2,H,theta1,theta2,V,sigma);

%
X1 = bridge.X1;
X2 = bridge.X2;

np = 51;

Fdirect = AxisymEvolverDirectForce(bridge,np,[],[]);
% Fderiv = AxisymEvolverDerivForce(bridge,np,[],[]);

Fdirect
% Fderiv






