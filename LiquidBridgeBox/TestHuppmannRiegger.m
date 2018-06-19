
clear all;

R = 1.0;
Vsph = SphereVolume(R);
sigma = 1.0;

%
% V = 0.1*Vsph;
% theta = deg2rad(36);
% V = 0.2*Vsph;
% theta = deg2rad(36);
% V = 0.1*Vsph;
% theta = deg2rad(90);
V = 0.2*Vsph;
theta = deg2rad(90);

%
dist = 0:0.1:1.2;
force = [];
for H = dist
	force(end+1) = HuppmannRiegger(R,H,theta,V,sigma);
end

if 1
	figure;
	plot(dist,force);
	axis([0,1,-3,3]);
end

