
addpath('../CubeChoppingBox');


D = 1.0;
R = D / 2;
Lx = 10.0;
Ly = 10.0;
xlo = 0; xhi = Lx;
ylo = 0; yhi = Ly;

dh = D / 30.0;
nx = round(Lx / dh);
ny = round(Ly / dh);
dx = Lx / nx;
dy = Ly / ny;

lambda = D / 1.5;
kappa = 1.0 / lambda;
kappa2 = kappa^2;

cellx = linspace(xlo+dx/2,xhi-dx/2,nx);
celly = linspace(ylo+dy/2,yhi-dy/2,ny);
nodex = linspace(xlo,xhi,nx+1);
nodey = linspace(ylo,yhi,ny+1);
[xcell,ycell] = ndgrid(cellx,celly);
[xumac,yumac] = ndgrid(nodex,celly);
[xvmac,yvmac] = ndgrid(cellx,nodey);

% setup particles
npart = 3;
partpos(1,:) = [ 2.0, 3.0 ];
partrad(1) = 1;
partdu(1) = 3.0;
partpos(2,:) = [ 5.0, 7.0 ];
partrad(2) = 2.0;
partdu(2) = 2.0;
partpos(3,:) = [ 8.0, 4.0 ];
partrad(3) = 1.5;
partdu(3) = 1.0;

% SDF and fraction
sdf = zeros(nx,ny); sdf(:) = 99999;
frac = zeros(nx,ny);
for ipart = 1:npart
    xcen = partpos(ipart,1);
    ycen = partpos(ipart,2);
    rcen = partrad(ipart);
    ducen = partdu(ipart);
    
    partsdf = CalcPartSDF(xcen,ycen,rcen, xcell,ycell);
    sdf = min(sdf,partsdf);
    
    partfrac = CalcPartFrac(xcen,ycen,rcen, xcell,ycell,dx,dy);
    frac = frac + partfrac;    
end

% tag for easy access different regions
tag = GhostIBTag(nx,ny,sdf);
tag_fluid = find(tag==1);
tag_nonfd = find(tag~=1);
tag_solid = find(tag==-1);
tag_ghost = find(tag==0);


% BC
bctype = [ 2, 2; 1, 1 ];
bcval = [ 0.0, 0.0; 0.0, 0.0 ];

% for matrix indexing
ind = reshape(1:nx*ny, nx,ny);

