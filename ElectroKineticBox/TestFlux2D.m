
clear all;

D = 1.0;
R = D / 2;
Lx = 10.0;
Ly = 10.0;
xlo = 0; xhi = Lx;
ylo = 0; yhi = Ly;


dh = D / 16.0;
nx = round(Lx / dh);
ny = round(Ly / dh);
dx = Lx / nx;
dy = Ly / ny;


xcell = linspace(xlo+dx/2,xhi-dx/2,nx);
ycell = linspace(ylo+dy/2,yhi-dy/2,ny);
[xcell,ycell] = ndgrid(xcell,ycell);

npart = 2;
partpos(1,:) = [ 3.0, 5.0 ];
partrad(1) = 1;
partdu(1) = -4.0;
partpos(2,:) = [ 7.0, 5.0 ];
partrad(2) = 1.0;
partdu(2) = -4.0;

sdf = zeros(nx,ny); sdf(:) = 99999;
phi = zeros(nx,ny);
dphi = zeros(nx,ny);
bphi = zeros(nx,ny);
for ipart = 1:npart
    xcen = partpos(ipart,1);
    ycen = partpos(ipart,2);
    partsdf = sqrt((xcell-xcen).^2 + (ycell-ycen).^2) - partrad(ipart);
    
    thickness = dh * 0.6;
    thd = tanh(-partsdf ./ thickness);
    partphi = 0.5 .* thd + 0.5;
    partdphi = -1/(2*thickness) .* (1.0 - thd.^2);
    partdphi = -partdphi;
    
    sdf = min(sdf,partsdf);
    phi = phi + partphi;
    dphi = dphi + partdphi;
    bphi = bphi + partdu(ipart).*partdphi;
end


bctype = [ 2, 2; 1, 1 ];
bcval = [ 0.0, 0.0; 0.0, 0.0 ];

tic;
[ ALap, rLap ] = MakeLap2Da(nx,ny,dx,dy, bctype,bcval);
toc;


rhs = -rLap + bphi(:);

sol = ALap \ rhs;
sol = reshape(sol,nx,ny);

figure;
% imagesc(xcell,ycell,sol);
contourf(xcell,ycell,sol);
% contourf(xcell,ycell,phi);
% contourf(xcell,ycell,dphi);
axis equal;
hold on;
contour(xcell,ycell,sdf,[0,0]);
hold off;


