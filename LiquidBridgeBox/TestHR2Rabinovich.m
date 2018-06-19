%%
%% (Rabinovich et al., 2005)
%%

clear;

% CA
theta = deg2rad(10);

% case 1
sigma = 27; % mN/m
R1 = 19; % um
R2 = 35; 
V = 0.2; % um^3
dist = 0:0.01:0.77;

% % case 2
% sigma = 24; % mN/m
% R1 = 19; % um
% R2 = 32.5; 
% V = 1.2; % um^3
% dist = 0:0.01:1.25;

% % case 3
% sigma = 28; % mN/m
% R1 = 19; % um
% R2 = 27.5; 
% V = 3.6; % um^3
% dist = 0:0.01:2.0;


data0 = [];
data1 = [];

if 1
    % Model 
	for H = dist
		data1(end+1) = BridgeForceHR2(R1*1e-6,R2*1e-6,H*1e-6,theta,theta*1e-18,V*1e-18,sigma);
	end
end

if 1
    % Optimal solution
	for H = 0:0.05:dist(end)
		np = 51;
		bridge = MakeBridge(R1,R2,H,theta,theta,V,sigma*1e-3);
		data0(end+1) = AxisymEvolverDirectForce(bridge,np,[],[]);
	end
end







