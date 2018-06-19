
%%
%% 
%%

clear all;

sigma = 1;

% R1 = 1.0;
% R2 = 2.0;
% theta1 = deg2rad(0);
% theta2 = deg2rad(50);
% V = 0.005;
% dist = linspace(0,0.1,21);

% R1 = 1.0;
% R2 = 2.0;
% theta1 = deg2rad(15);
% theta2 = deg2rad(60);
% V = 0.02;
% dist = linspace(0,0.2,21);

% R1 = 1.0;
% R2 = 2.0;
% theta1 = deg2rad(15);
% theta2 = deg2rad(60);
% V = 0.04;
% dist = linspace(0,0.2,21);
R1 = 1.0;
R2 = 3.0;
theta1 = deg2rad(15);
theta2 = deg2rad(60);
V = 0.04;
dist = linspace(0,0.2,21);

% R1 = 1.0;
% R2 = 2.0;
% theta1 = deg2rad(120);
% theta2 = deg2rad(150);
% V = 0.04;
% dist = linspace(0,0.35,8);
% R1 = 2.0;
% R2 = 1.0;
% theta1 = deg2rad(120);
% theta2 = deg2rad(150);
% V = 0.04;
% dist = linspace(0,0.35,8);


%% bi-concave
% R1 = 1.0;
% R2 = 2.0;
% theta1 = deg2rad(30);
% theta2 = deg2rad(60);
% V = 0.04;
% dist = linspace(0,0.35,8);

%% bi-convex
% R1 = 1.0;
% R2 = 2.0;
% theta1 = deg2rad(120);
% theta2 = deg2rad(150);
% V = 0.04;
% dist = linspace(0,0.35,8);

data0 = [];
data1 = [];
data2 = [];
data3 = [];

if 1
	for H = dist
		Rder = DerjaguinRadius(R1,R2);
		tder = DerjaguinAngle(theta1,theta2);
		data1(end+1) = BridgeForceRabinovich(Rder,H,tder,V,sigma);
		data2(end+1) = BridgeForceHuppmannRiegger(Rder,H,tder,V,sigma);
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

data = [data0',data1',data2',data3'];

