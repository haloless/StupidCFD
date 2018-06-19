
clear;

%% TOYOTA, particle-wall bridge force


sigma = 1.0;

R1 = 1.0;
R2 = -1.0; % wall

theta1 = deg2rad(60);
theta2 = deg2rad(100);

% the average CA
thetam = acos(0.5*(cos(theta1)+cos(theta2)));

% 40% liquid volume
Vliq = 4/3*pi * R1^3 * 0.4;


% ncoord = 1;
% ncoord = 2;
% ncoord = 4;
% ncoord = 6;
% ncoord = 8;
ncoord = 10;

V = Vliq / ncoord;


data = [];
for H = 0:0.01:0.5
    f1 = BridgeForceHR2(R1,R2,H,theta1,theta2,V,sigma);
    f2 = BridgeForceHR2(R1,R2,H,thetam,thetam,V,sigma);
    data(end+1,:) = [H, f1,f2];
end









