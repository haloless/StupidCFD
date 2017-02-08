function [sdf] = MakeSDF(shape)


lx = shape.a * 6;
ly = shape.b * 6;
xmin = shape.xc - lx/2;
xmax = shape.xc + lx/2;
ymin = shape.yc - ly/2;
ymax = shape.yc + ly/2;

dh = min(shape.a,shape.b) / 20;
[xg,yg] = ndgrid(xmin:dh:xmax,ymin:dh:ymax);
nx = size(xg,1);
ny = size(xg,2);

phig = zeros(nx,ny);
for i = 1:nx
for j = 1:ny
	phig(i,j) = ShapePotential(shape,xg(i,j),yg(i,j));
end
end

% reinitialize to get SDF
band = lx * 0.2;
phig = ImplicitFuncReinit(phig,nx,ny,dh,dh, band);


sdf = struct();
sdf.xmin = xmin;
sdf.xmax = xmax;
sdf.ymin = ymin;
sdf.ymax = ymax;
sdf.dh = dh;
sdf.nx = nx;
sdf.ny = ny;
sdf.xg = xg;
sdf.yg = yg;
sdf.phig = phig;



return
end

