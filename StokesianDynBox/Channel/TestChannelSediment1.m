
clear all;

addpath('../');

CommonGlobals;
CommonInit;

ChannelMobGlobals;
ChannelMobLoadAll;

a = 1.0;
Fext = [ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]';
ewall = [ 0.0, 0.0, 1.0 ]';

Hs = [10, 15, 20, 30, 50, 75, 100];
% Hs = [10];

ts = [];
ts = [ts, linspace(1.0e-5,9.0e-5,9)];
ts = [ts, linspace(1.0e-4,9.0e-4,9)];
ts = [ts, linspace(1.0e-3,9.0e-3,9)];
ts = [ts, linspace(1.0e-2,9.0e-2,9)];
ts = [ts, linspace(1.0e-1,5.0e-1,5)];

vel = [];
ang = [];

for H = Hs
	ah = a / H;
	
	tmp = [];
	tmp2 = [];
	for t = ts
		theta = t + ah;
		
		muf = PartChanMobUF(theta,ah,ewall);
		mul = PartChanMobUL(theta,ah,ewall);
		mol = PartChanMobOL(theta,ah,ewall);
		
		mat = [ muf, mul; mul',mol ];
		
		sol = mat * Fext;
		tmp(end+1) = sol(1);
		tmp2(end+1) = sol(5);
	end
	vel(:,end+1) = tmp';
	ang(:,end+1) = tmp2';
end

figure;
loglog(ts',vel); 
axis([10^-5 0.5 10^-1 10^0]);

figure;
semilogx(ts',ang); 
axis([10^-6 0.5 -0.01 0.05]);

