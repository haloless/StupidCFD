
clear all;

% gap size
hs = [];
hs = [ hs, linspace(0.00001,0.00009,9)];
hs = [ hs, linspace(0.0001,0.0009,9)];
hs = [ hs, linspace(0.001,0.009,9)];
hs = [ hs, linspace(0.01,0.09,9)];
hs = [ hs, linspace(0.1,0.9,9)];
hs = [ hs, linspace(1.0,3.0,3)];

% parallel force
f1 = [1,0,0, 0,0,0]';
% parallel velocity
ux1 = [];
oy1 = [];

% vertical force
f2 = [0,0,1, 0,0,0]';
uz2 = [];

for h = hs
	% center position
	z = h + 1;
	
	res = WallResistanceMatrix(z);
	
	% parallel
	vel = res \ f1;
	ux1(end+1) = vel(1);
	oy1(end+1) = vel(5);
	
	% vertical
	vel = res \ f2;
	uz2(end+1) = vel(3);
end

figure;
loglog(hs./100, ux1);
axis([10^-5, 0.5, 10^-1, 10^0]);

figure;
semilogx(hs./100, oy1);
axis([10^-5, 0.5, -0.01,0.05]);

figure;
loglog(hs./100, uz2);
% axis([10^-3, 0.5, 10^-4, 10^-1]);




