%%
%% (Lambert, 2008)
%%

clear all;

sigma = 0.0208;
R1 = 2.0;
R2 = 3.95;
theta1 = deg2rad(0);
theta2 = deg2rad(14.3);

% case 1, small volume
% V = 0.065;
% dist = linspace(0,0.32,17);
% case 2, mediate volume
V = 1.4;
dist = linspace(0,0.60,13);

data0 = [];
data1 = [];
data2 = [];
data3 = [];

if 1
	for H = dist
		R = DerjaguinRadius(R1,R2);
		theta = DerjaguinAngle(theta1,theta2);
		data1(end+1) = BridgeForceRabinovich(R,H,theta,V,sigma);
		data2(end+1) = BridgeForceHuppmannRiegger(R,H,theta,V,sigma);
		data3(end+1) = BridgeForceHR2(R1,R2,H,theta1,theta2,V,sigma);
	end
end

if 1
	for H = dist
		np = 51;
		bridge = MakeBridge(R1,R2,H,theta1,theta2,V,sigma);
		data0(end+1) = AxisymEvolverDirectForce(bridge,np,[],[]);
	end
end

data = [dist', data0',data1',data2',data3'];






