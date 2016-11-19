
clear all;

D = 1.0;
R = D / 2;
Lx = 10.0;
Ly = 10.0;
xlo = 0; xhi = Lx;
ylo = 0; yhi = Ly;


dh = D / 10.0;
nx = round(Lx / dh);
ny = round(Ly / dh);
dx = Lx / nx;
dy = Ly / ny;


xcell = linspace(xlo+dx/2,xhi-dx/2,nx);
ycell = linspace(ylo+dy/2,yhi-dy/2,ny);
[xcell,ycell] = ndgrid(xcell,ycell);

npart = 3;
partpos(1,:) = [ 2.0, 3.0 ];
partrad(1) = 1;
partdu(1) = -3.0;
partpos(2,:) = [ 5.0, 7.0 ];
partrad(2) = 2.0;
partdu(2) = -2.0;
partpos(3,:) = [ 8.0, 4.0 ];
partrad(3) = 1.5;
partdu(3) = -1.0;

sdf = zeros(nx,ny); sdf(:) = 99999;
owner = zeros(nx,ny);
for ipart = 1:npart
    xcen = partpos(ipart,1);
    ycen = partpos(ipart,2);
    partsdf = sqrt((xcell-xcen).^2 + (ycell-ycen).^2) - partrad(ipart);
    
    sdf = min(sdf,partsdf);
    owner(partsdf<=0) = ipart;
end

tag = zeros(nx,ny);
tag(sdf>0) = 1;
tag(sdf<=0) = -1;
for j = 2:ny-1
for i = 2:nx-1
    if sdf(i,j) <= 0
        if sdf(i-1,j)>0 | sdf(i+1,j)>0 | sdf(i,j-1)>0 | sdf(i,j+1)>0
            tag(i,j) = 0;
        end
    end
end
end
tag_fluid = find(tag==1);
tag_solid = find(tag==-1);
tag_ghost = find(tag==0);


bctype = [ 2, 2; 1, 1 ];
bcval = [ 0.0, 0.0; 0.0, 0.0 ];

ind = reshape(1:nx*ny, nx,ny);

disp(['Build Laplacian']);
tic;
[ ALap, rLap ] = MakeLap2Da(nx,ny,dx,dy, bctype,bcval);
toc;

disp(['Build Interpolation'])
Aint = spdiags(ones(nx*ny,1), 0, nx*ny,nx*ny);
bint = zeros(nx*ny,1);
tic;
for j = 1:ny
for i = 1:nx
if tag(i,j) == 0
    dist = abs(sdf(i,j));
    nvecx = (sdf(i+1,j)-sdf(i-1,j)) / (dx*2);
    nvecy = (sdf(i,j+1)-sdf(i,j-1)) / (dy*2);
    nnorm = sqrt(nvecx^2 + nvecy^2);
    nvecx = nvecx / nnorm;
    nvecy = nvecy / nnorm;
    
    % this point
    xx = xcell(i,j);
    yy = ycell(i,j);
    
    % on boundary
    xb = xx + dist*nvecx;
    yb = yy + dist*nvecy;
    
    % control point
    len = dh * 1.5;
    l1 = len * 1;
    x1 = xb + l1*nvecx;
    y1 = yb + l1*nvecy;
    l2 = len * 2;
    x2 = xb + l2*nvecx;
    y2 = yb + l2*nvecy;
    
    c0 = -(l1+l2) / (dist+l1) / (dist+l2);
    c1 = (dist-l2) / (dist+l1) / (l1-l2);
    c2 = (dist-l1) / (dist+l2) / (l2-l1);
    
    dub = partdu(owner(i,j));
    
    idx = ind(i,j);
    Aint(idx,idx) = 0;
    bint(idx) = dub / c0;
    
    [ is, js, ws ] = CellInterpCoef2D(x1,y1, nx,ny,dx,dy,xlo,ylo);
    for k = 1:4
        idk = ind(is(k),js(k));
        if tag(idk) ~= 1
            error('bad interp');
        end
        Aint(idx,idk) = Aint(idx,idk) + c1*ws(k)/(-c0);
    end
    [ is, js, ws ] = CellInterpCoef2D(x2,y2, nx,ny,dx,dy,xlo,ylo);
    for k = 1:4
        idk = ind(is(k),js(k));
        if tag(idk) ~= 1
            error('bad interp');
        end
        Aint(idx,idk) = Aint(idx,idk) + c2*ws(k)/(-c0);
    end
end
end
end
toc;



rhs = -rLap;
sol = zeros(nx*ny,1);
% diagonal of ALap, needed by relaxation
dLap = spdiags(ALap, 0);

% initial
sol = Aint*sol + bint;
res = rhs - ALap*sol;
res(tag_ghost) = 0;
rnorm0 = norm(res);
disp(['|res0|=',num2str(rnorm0)]);

if (1)
% use simple relaxation
for iter = 1:99999
    sold = sol;
    sol = (rhs - ALap*sol) ./ dLap + sol;
    
    omega = 1.0;
    sol = omega.*sol + (1-omega).*sold;
    
    
    sol = Aint*sol + bint;
    
    res = rhs - ALap*sol;
    res(tag_ghost) = 0;
    rnorm = norm(res);
    
    if mod(iter,100) == 0
        disp(['|res|=',num2str(rnorm)]);
    end
    if rnorm <= rnorm0*1e-8 | rnorm<=1e-8
        disp(['converged']);
        disp(['|res|=',num2str(rnorm)]);
        break;
    end
end
end

sol = reshape(sol,nx,ny);

sol(tag_solid) = nan;
% sol(tag_ghost) = nan;

figure;
% imagesc(xcell,ycell,sol);
contourf(xcell,ycell,sol);
% contourf(xcell,ycell,phi);
% contourf(xcell,ycell,dphi);
% imagesc([xlo+dx/2,xhi-dx/2],[ylo+dy/2,yhi-dy/2],owner');
% imagesc([xlo+dx/2,xhi-dx/2],[ylo+dy/2,yhi-dy/2],tag');
axis equal;
hold on;
contour(xcell,ycell,sdf,[0,0]);
hold off;

% analytical flux from internal boundary
% = (particle circumference) * (flux density)
flux_ana = sum(2*pi*partrad.*partdu)
if (1)
    % check numerical flux through external wall 
    % they should be equal
    flux_wall = 0;
    if bctype(2,1) == 1
        flux_wall = flux_wall + sum(bcval(2,1)-sol(:,1))/(dy/2)*dx;
    end
    if bctype(2,2) == 1
        flux_wall = flux_wall + sum(bcval(2,2)-sol(:,ny))/(dy/2)*dx;
    end
    flux_wall
end
if (1)
    % check numerical flux through internal particles
    % they should be equal to the numerical flux through external wall
    flux_part = 0;
    for j = 1:ny
    for i = 1:nx
    if tag(i,j) == 0
        if tag(i+1,j) == 1
            ff = (sol(i+1,j)-sol(i,j)) / dx * dy;
            flux_part = flux_part + ff;
        end
        if tag(i-1,j) == 1
            ff = (sol(i-1,j)-sol(i,j)) / dx * dy;
            flux_part = flux_part + ff;
        end
        if tag(i,j+1) == 1
            ff = (sol(i,j+1)-sol(i,j)) / dy * dx;
            flux_part = flux_part + ff;
        end
        if tag(i,j-1) == 1
            ff = (sol(i,j-1)-sol(i,j)) / dy * dx;
            flux_part = flux_part + ff;
        end
    end
    end
    end
    flux_part
end












