
clear all;

sigma = 1.0;
R = 1.0;

% V = SphereVolume(R)*0.1;
% dist = linspace(0,0.5,11);
V = SphereVolume(R)*0.2;
dist = linspace(0,1.0,11);

% theta = deg2rad(15);
% theta = deg2rad(45);
% theta = deg2rad(90);
theta = deg2rad(120);


data0 = [];
data0a = [];
data1 = [];
data2 = [];

if 1
	for H = dist
		% data1(end+1) = BridgeForceRabinovich(R,H,theta,V,sigma);
		% data1(end+1) = BridgeForceHuppmannRiegger(R,H,theta,V,sigma);
		data2(end+1) = BridgeForceHR2(R,R,H,theta,theta,V,sigma);
	end
end

if 0
	for H = dist
		np = 51;
		bridge = MakeBridge(R,R,H,theta,theta,V,sigma);
		[data0(end+1),~,~, data0a(end+1),~,~] = AxisymEvolverDirectForce(bridge,np,[],[]);
	end
end

data = [data0',data1',data2'];






