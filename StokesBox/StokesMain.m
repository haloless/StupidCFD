
clear all;

%
ImUnit = 0.0 + 1.0i;

Re = 10.0;

Lx = 2*pi;
Ly = 2*pi;

% nn = 32;
nn = 64;
nx = nn;
ny = nn;
dx = Lx / nx;
dy = Ly / ny;

xs = linspace(0,Lx-dx,nx)';
ys = linspace(0,Ly-dy,ny)';
[xs,ys] = ndgrid(xs,ys);

kx = zeros(nx,ny);
ky = zeros(nx,ny);
for j = 1:ny
	kx(:,j) = [0:nx/2-1, -nx/2:-1]';
end
for i = 1:nx
	% ky(i,:) = [0:-1:-ny/2+1, ny/2:-1:1];
	ky(i,:) = [0:ny/2-1, -ny/2:-1];
end

k2 = kx.^2 + ky.^2;
% k2(1,1) = 1;

fun_gx = @(x,y) exp(-4*((x-pi).^2+(y-pi).^2)) .* (2+tanh(y-pi));
gx = fun_gx(xs,ys);
gy = zeros(nx,ny);


figure;
quiver(xs,ys,gx,gy);
axis equal;
axis([0 Lx 0 Ly]);

gxhat = fft2(gx);
gyhat = fft2(gy);

if (0)
	% manually inverse transform
	% test configuration of wave number
	gx2r = zeros(nx,ny);
	gx2m = zeros(nx,ny);
	for j = 1:ny
	for i = 1:nx
		gij = 0;
		for jj = 1:ny
		for ii = 1:nx
			kbyx = kx(ii,jj)*xs(i,j) + ky(ii,jj)*ys(i,j);
			gij = gij + gxhat(ii,jj)*exp(ImUnit*kbyx);
		end
		end
		gij = gij / (nx*ny);
		gx2r(i,j) = real(gij);
		gx2m(i,j) = imag(gij);
	end
	end
	grerr = max(abs(gx(:)-gx2r(:)))
	gierr = max(abs(gx2m(:)))
end

% in wave-number space
uhat = zeros(nx,ny);
vhat = zeros(nx,ny);
phat = zeros(nx,ny);
% global Grad(p) to balance force
pgx = 0;
pgy = 0;

for jj = 1:ny
for ii = 1:nx
	ki = kx(ii,jj);
	kj = ky(ii,jj);
	kk = k2(ii,jj);
	
	gi = gxhat(ii,jj);
	gj = gyhat(ii,jj);
	
	if kk == 0
		uhat(ii,jj) = 0;
		vhat(ii,jj) = 0;
		phat(ii,jj) = 0;
		pgx = gi;
		pgy = gj;
	else
		uhat(ii,jj) = 1/kk / (kk/Re) * (kj*kj*gi - ki*kj*gj);
		vhat(ii,jj) = 1/kk / (kk/Re) * (-ki*kj*gi + ki*ki*gj);
		phat(ii,jj) = 1/kk * (-ImUnit) * (ki*gi + kj*gj);
	end
end
end

us = real(ifft2(uhat));
vs = real(ifft2(vhat));
ps = real(ifft2(phat));


xper = linspace(0,Lx,nx+1);
yper = linspace(0,Ly,ny+1);
[xper,yper] = ndgrid(xper,yper);
pper = zeros(nx+1,ny+1);
pper(1:nx,1:ny) = ps;
pper(nx+1,:) = pper(1,:);
pper(:,ny+1) = pper(:,1);

figure;
% contourf(xs,ys,ps);
contourf(xper,yper,pper,24);
hold on;
quiver(xs,ys,us,vs);
hold off;
axis equal;
axis([0 Lx 0 Ly]);

figure; surf(xper,yper,pper)


