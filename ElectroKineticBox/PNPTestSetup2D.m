
% addpath('../CubeChoppingBox');


%
% domain
%
D = 1.0;
R = D / 2;
Lx = 5.0;
Ly = 16.0;
% Lx = 20.0;
% Ly = 20.0;
xlo = 0; xhi = Lx;
ylo = 0; yhi = Ly;

%
% setup grid
%
dh = D / 10.0;
nx = round(Lx / dh);
ny = round(Ly / dh);
dx = Lx / nx;
dy = Ly / ny;

cellx = linspace(xlo+dx/2,xhi-dx/2,nx);
celly = linspace(ylo+dy/2,yhi-dy/2,ny);
nodex = linspace(xlo,xhi,nx+1);
nodey = linspace(ylo,yhi,ny+1);
[xcell,ycell] = ndgrid(cellx,celly);
[xumac,yumac] = ndgrid(nodex,celly);
[xvmac,yvmac] = ndgrid(cellx,nodey);

%
% setup problem
%

%
lambda = 2;
clambda2 = 0.5*lambda^2;

if (1)
% imposed electric field
E0 = [ 0.3865, 0.0 ];

% setup particles
npart = 1;
partpos(1,:) = [ 2.5, 1.0 ];
% partpos(1,:) = [ 2.5, 10.0 ];
partrad(1) = 0.5;

% external BC
bctype_ext = [ 2, 2; 0, 0 ];
bcval_ext = [ -E0(1)*xlo, -E0(1)*xhi; 0.0, 0.0 ];
partbc_ext(1) = 0;
partval_ext(1) = 0.0;
% internal BC
bctype_int = [ 2, 2; 1, 1 ];
bcval_int = [ 0.0, 0.0; -1.0, 0.0 ];
partbc_int(1) = 1;
partval_int(1) = -1.0;
end

if (0)
% imposed electric field
E0 = [ 0.0, 0.0 ];

% setup particles
npart = 1;
partpos(1,:) = [ 10.0, 10.0 ];
partrad(1) = 0.5;

% external BC
bctype_ext = [ 1, 1; 1, 1 ];
bcval_ext = [ 0.0, 0.0; 0.0, 0.0 ];
partbc_ext(1) = 0;
partval_ext(1) = 0.0;
% internal BC
bctype_int = [ 1, 1; 1, 1 ];
bcval_int = [ 0.0, 0.0; 0.0, 0.0 ];
partbc_int(1) = 1;
partval_int(1) = -1.0;
end

% species
nspec = 2;
specz(1) = 1;
specz(2) = -1;
specdens = zeros(nx,ny,nspec);
specdens(:,:,1) = 1.0;
specdens(:,:,2) = 1.0;

% SDF and fraction
sdf = zeros(nx,ny); sdf(:) = 99999;
% frac = zeros(nx,ny);
owner = zeros(nx,ny);
for ipart = 1:npart
    xcen = partpos(ipart,1);
    ycen = partpos(ipart,2);
    rcen = partrad(ipart);
    
    partsdf = CalcPartSDF(xcen,ycen,rcen, xcell,ycell);
    sdf = min(sdf,partsdf);
    
    % partfrac = CalcPartFrac(xcen,ycen,rcen, xcell,ycell,dx,dy);
    % frac = frac + partfrac;
    
    owner(partsdf<=0) = ipart;
end

% tag for easy access different regions
tag = GhostIBTag(nx,ny,sdf);
tag_fluid = find(tag==1);
tag_nonfd = find(tag~=1);
tag_solid = find(tag==-1);
tag_ghost = find(tag==0);



% for matrix indexing
matind = reshape(1:nx*ny, nx,ny);

