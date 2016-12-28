
clear all;

mu = 1.0/(6*pi);

lx = pi*2;
ly = pi*2;
lz = pi*2;

f0 = -1.0;
alpha = 3.0;

nmax = 10;
kmax = 10;

dh = 0.001;

dist = [0.01:0.01:0.09, 0.1:0.1:3.0];

ux = [];
for x = dist
	y = 0;
	z = 0;
	
	s1 = EwaldS1(x,y,z, lx,ly,lz, alpha, nmax,kmax);
	
	s2 = EwaldS2(x,y,z, lx,ly,lz, alpha, nmax,kmax);
	s2p = EwaldS2(x+dh,y,z, lx,ly,lz, alpha, nmax,kmax);
	s2m = EwaldS2(x-dh,y,z, lx,ly,lz, alpha, nmax,kmax);
	ds2 = (s2p-2*s2+s2m) / dh^2;
	
	u = -f0/(4*pi*mu) * (s1-ds2);
	ux(end+1) = u;
end

uy = [];
for y = dist
	x = 0;
	z = 0;
	
	s1 = EwaldS1(x,y,z, lx,ly,lz, alpha, nmax,kmax);
	
	s2 = EwaldS2(x,y,z, lx,ly,lz, alpha, nmax,kmax);
	s2p = EwaldS2(x+dh,y,z, lx,ly,lz, alpha, nmax,kmax);
	s2m = EwaldS2(x-dh,y,z, lx,ly,lz, alpha, nmax,kmax);
	ds2 = (s2p-2*s2+s2m) / dh^2;
	
	u = -f0/(4*pi*mu) * (s1-ds2);
	uy(end+1) = u;
end



