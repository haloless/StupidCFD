%%
%% (Willet et al., 2000)
%%

clear all;

% fixed values
sigma = 0.0204;
theta = deg2rad(0);

%% R1/R2 = 1:1 
% R1 = 2.381;
% R2 = 2.381;

% V = 0.0136;
% dist = linspace(0,0.24,13);
% V = 0.0313;
% dist = linspace(0,0.32,17);
% V = 0.0742;
% dist = linspace(0,0.42,22);

%% R1/R2 = 1.5:1
% R1 = 2.381;
% R2 = 1.588;

% V = 0.0096;
% dist = linspace(0,0.22,12);
% V = 0.0132;
% dist = linspace(0,0.24,13);
% V = 0.0247;
% dist = linspace(0,0.30,16);
% V = 0.0593;
% dist = linspace(0,0.40,21);

%% R1/R2 = 2:1
% R1 = 2.381;
% R2 = 1.191;

% V = 0.0253;
% dist = linspace(0,0.30,16);
% V = 0.0618;
% dist = linspace(0,0.40,21);
% V = 0.1278;
% dist = linspace(0,0.50,26);


%% wall
R1 = 2.381;
% R2 = -1;
R2 = 1000; % large right sphere
% V = 0.1624;
% dist = linspace(0,0.5,21);
V = 0.2802;
dist = linspace(0,0.6,25);


data0 = [];
data1 = [];
data2 = [];
data3 = [];

if 1
	for H = dist
		R = DerjaguinRadius(R1,R2);
		% data1(end+1) = BridgeForceRabinovich(R,H,theta,V,sigma);
		% data2(end+1) = BridgeForceHuppmannRiegger(R,H,theta,V,sigma);
		data3(end+1) = BridgeForceHR2(R1,R2,H,theta,theta,V,sigma);
	end
end

if 0
	for H = dist
		% np = 81;
		np = 51;
		bridge = MakeBridge(R1,R2,H,theta,theta,V,sigma);
		data0(end+1) = AxisymEvolverDirectForce(bridge,np,[],[]);
	end
end

data = [dist', data0',data1',data2',data3'];






