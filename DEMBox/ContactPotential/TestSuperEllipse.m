
clear all;

nx = 33;
ny = 33;
nz = 33;
xs = linspace(-3,3,nx);
ys = linspace(-3,3,ny);
zs = linspace(-3,3,nz);

[xs,ys,zs] = ndgrid(xs,ys,zs);
% [xs,ys,zs] = meshgrid(xs,ys,zs);

a = 1.0;
b = 1.0;
c = 2.0;
alpha = 2.0;
beta = 2.0;
gamma = 16.0;

phis = (xs/a).^alpha + (ys/b).^beta + (zs/c).^gamma - 1;

% phis = zeros(nx,ny,nz);
% for i = 1:nx
% for j = 1:ny
% for k = 1:nz
	% phis(i,j,k) = (xs(i)/a)^alpha + (ys(j)/b)^beta + (zs(k)/c)^gamma - 1;
% end
% end
% end

figure;
% isosurface(xs,ys,zs,phis, 0);
isosurface(phis, 0);
axis equal;



