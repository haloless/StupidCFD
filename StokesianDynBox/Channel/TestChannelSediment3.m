
clear all;

CommonGlobals;
CommonInit;

ChannelMobGlobals;
ChannelMobLoadAll;

a = 1.0;
Fext = [ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]';
ewall = [ 0.0, 0.0, 1.0 ]';

% Hs = [10, 15, 20, 30, 50, 75, 100];
Hs = [30, 50, 75, 100];
% Hs = [10];

% ts = [];
% ts = [ts, linspace(1.0e-5,9.0e-5,9)];
% ts = [ts, linspace(1.0e-4,9.0e-4,9)];
% ts = [ts, linspace(1.0e-3,9.0e-3,9)];
% ts = [ts, linspace(1.0e-2,9.0e-2,9)];
% ts = [ts, linspace(1.0e-1,5.0e-1,5)];

ts = 0.05:0.01:0.5;

vel = [];

for H = Hs
	ah = a / H;
	
	tmp = [];
	for t = ts
		theta = t;
		
		muf = PartChanMobUF(theta,ah,ewall);
		mul = PartChanMobUL(theta,ah,ewall);
		mol = PartChanMobOL(theta,ah,ewall);
		
		mat = [ muf, mul; mul',mol ];
		
		sol = mat * Fext;
		tmp(end+1) = sol(1);
	end
	vel(:,end+1) = tmp';
end

figure;
plot(ts',vel,'x-'); 
axis([0,1.0, 0.5,1]);
% axis([10^-5 0.5 10^-1 10^0]);





