
clear all;

PBTestSetup2D;

owner = zeros(nx,ny);
for ipart = 1:npart
    xcen = partpos(ipart,1);
    ycen = partpos(ipart,2);
    partsdf = sqrt((xcell-xcen).^2 + (ycell-ycen).^2) - partrad(ipart);
    
    owner(partsdf<=0) = ipart;
end




disp(['Build Laplacian']);
tic;
[ ALap, rLap ] = MakeLap2Da(nx,ny,dx,dy, bctype,bcval);
toc;

bint = zeros(nx*ny,1);
if (1)
disp(['Build Interpolation'])
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
    ALap(idx,:) = 0;
    ALap(idx,idx) = c0;
    bint(idx) = dub;
    
    [ is, js, ws ] = CellInterpCoef2D(x1,y1, nx,ny,dx,dy,xlo,ylo);
    for k = 1:4
        idk = ind(is(k),js(k));
        if tag(idk) ~= 1
            error('bad interp');
        end
        ALap(idx,idk) = ALap(idx,idk) + c1*ws(k);
    end
    [ is, js, ws ] = CellInterpCoef2D(x2,y2, nx,ny,dx,dy,xlo,ylo);
    for k = 1:4
        idk = ind(is(k),js(k));
        if tag(idk) ~= 1
            error('bad interp');
        end
        ALap(idx,idk) = ALap(idx,idk) + c2*ws(k);
    end
end
end
end
toc;
end

sol = zeros(nx*ny,1);
% initial residual
bpb = kappa2.*sinh(sol);
bpb(tag_nonfd) = 0;
res = ALap*sol + rLap - bint - bpb;
res(tag_solid) = 0;
rnorm0 = norm(res);
disp(['|res0|=',num2str(rnorm0)]);
for cycle = 1:100
    
    % poission-boltzmann source term
    % this term only applies to true fluid cells
    % so the non-fluid parts are set to zero
    Apb = kappa2 .* cosh(sol);
    Apb(tag_nonfd) = 0;
    bpb = kappa2.*sinh(sol) - kappa2.*cosh(sol).*sol;
    bpb(tag_nonfd) = 0;
    Apb = spdiags(Apb,0, nx*ny,nx*ny);
    
    
    A = ALap - Apb;
    rhs = -rLap + bint + bpb;
    
    sol = A \ rhs;
    
    % check residual for nonlinear eq
    bpb = kappa2.*sinh(sol);
    bpb(tag_nonfd) = 0;
    res = ALap*sol + rLap - bint - bpb;
    res(tag_solid) = 0;
    rnorm = norm(res);
    disp(['cycle=',int2str(cycle), ', |res|=',num2str(rnorm)]);
    if rnorm <= rnorm0*1e-8
        disp(['converged']);
        break;
    end
end

sol = reshape(sol,nx,ny);

sol(tag_solid) = nan;
% sol(tag_ghost) = nan;

figure;
PBTestPlotSol;

PBTestPlotIon;

PBTestCheckSum;













