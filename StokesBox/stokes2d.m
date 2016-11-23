
clear all;

twopi = pi * 2;
ImUnit = 1i;

rho = 1.0;
nu = 1.0;

Lx = pi * 2;
Ly = pi * 2;

nx = 64;
ny = 64;
dx = Lx / nx;
dy = Ly / ny;

xs = linspace(0,Lx-dx,nx);
ys = linspace(0,Ly-dy,ny);
[xs,ys] = ndgrid(xs,ys);

xper = linspace(0,Lx,nx+1);
yper = linspace(0,Ly,ny+1);
[xper,yper] = ndgrid(xper,yper);

kx = zeros(nx,ny);
ky = zeros(nx,ny);
for j = 1:ny
    kx(:,j) = [0:nx/2-1, -nx/2:-1];
end
for i = 1:nx
    ky(i,:) = [0:ny/2-1, -ny/2:-1];
end
k2 = kx.^2 + ky.^2;

fs = zeros(nx,ny);
gs = zeros(nx,ny);
% fs(:,:) = sin(xs) .* cos(xs) + 2*nu*cos(xs).*sin(ys);
% gs(:,:) = sin(ys) .* cos(ys) - 2*nu*sin(xs).*cos(ys);
fs(:,round(ny*1/4)) = -1;
fs(:,round(ny*3/4)) = 1;

if (1)
    figure;
    quiver(xs,ys,fs,gs);
    axis equal;
    axis([0 Lx 0 Ly]);
end

fh = fft2(fs);
gh = fft2(gs);

%
uh = zeros(nx,ny);
vh = zeros(nx,ny);
ph = zeros(nx,ny);

for j = 1:ny
for i = 1:nx
    ki = kx(i,j);
    kj = ky(i,j);
    kn = k2(i,j);
    fhat = fh(i,j);
    ghat = gh(i,j);
    
    if kn == 0
        uh(i,j) = 0;
        vh(i,j) = 0;
        ph(i,j) = 0;
    else
        uh(i,j) = 1.0/kn / (kn*nu) * (kj*kj*fhat - ki*kj*ghat);
        vh(i,j) = 1.0/kn / (kn*nu) * (-ki*kj*fhat + ki*ki*ghat);
        ph(i,j) = 1.0/kn * (-ImUnit*ki*fhat - ImUnit*kj*ghat);
    end
end
end

%
us = real(ifft2(uh));
vs = real(ifft2(vh));
ps = real(ifft2(ph));

pper = zeros(nx+1,ny+1);
pper(1:nx,1:ny) = ps;
pper(:,ny+1) = pper(:,1);
pper(nx+1,:) = pper(1,:);

if 1
    figure;
    contourf(xper,yper,pper);
    hold on;
    quiver(xs,ys,us,vs);
    hold off;
    axis equal;
    axis([0 Lx 0 Ly]);
end

uana = cos(xs).*sin(ys);
vana = -sin(xs).*cos(ys);
pana = -0.25 * (cos(2*xs) + cos(2*ys));
