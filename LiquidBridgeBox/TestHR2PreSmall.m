
clear all;

sigma = 1.0;
R = 1.0;

% theta = deg2rad(15);
% theta = deg2rad(30);
% theta = deg2rad(60);
theta = deg2rad(90);

V = SphereVolume(R)*0.005;

dist = linspace(0,0.2,21);
data0 = [];
data1 = [];
data2 = [];
data3 = [];

if 1
	for H = dist
		% data1(end+1) = BridgeForceRabinovich(R,H,theta,V,sigma);
		% data2(end+1) = BridgeForceHuppmannRiegger(R,H,theta,V,sigma);
		data3(end+1) = BridgeForceHR2(R,R,H,theta,theta,V,sigma);
	end
end

if 0
	for H = dist
		np = 51;
		bridge = MakeBridge(R,R,H,theta,theta,V,sigma);
		data0(end+1) = AxisymEvolverDirectForce(bridge,np,[],[]);
	end
end

data = [data0',data1',data2',data3'];






