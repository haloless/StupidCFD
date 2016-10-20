
clear all;

PBTestSetup2D;


dfrac = zeros(nx,ny);
bsurf = zeros(nx,ny);
for ipart = 1:npart
    xcen = partpos(ipart,1);
    ycen = partpos(ipart,2);
    rcen = partrad(ipart);
    ducen = partdu(ipart);
    
    partfu = CalcPartFrac(xcen,ycen,rcen, xumac,yumac,dx,dy);
    partfv = CalcPartFrac(xcen,ycen,rcen, xvmac,yvmac,dx,dy);
    partdf = sqrt(((partfu(2:nx+1,:)-partfu(1:nx,:))./dx).^2 + ((partfv(:,2:ny+1)-partfv(:,1:ny))./dy).^2);
    % partdf = (1.0-partfrac) .* partfrac;
    
    dfrac = dfrac + partdf;
    bsurf = bsurf + partdf .* ducen;
end

if (0)
    figure;
    % contour(xcell,ycell,frac,[0.5,0.5]);
    % contourf(xcell,ycell,frac);
    contourf(xcell,ycell,bsurf);
    % contour(xcell,ycell,frac,[0.5,0.5]);
    axis equal;
    hold on;
    contour(xcell,ycell,frac,[0.5,0.5]);
    contour(xcell,ycell,sdf,[0.0,0.0]);
hold off;
end


disp(['Build Laplacian']);
tic;
[ ALap, rLap ] = MakeLap2Da(nx,ny,dx,dy, bctype,bcval);
toc;

% initial residual
sol = zeros(nx*ny,1);
bsurf = reshape(bsurf,nx*ny,1);
bfrac = reshape(1.0-frac, nx*ny,1);
bfrac(:) = 1; 
% bfrac(frac>0) = 0; 
bfrac(frac==1) = 0;
bfrac(dfrac>0) = 0;
bpb = kappa2.*sinh(sol).*bfrac;
res = ALap*sol + rLap - bsurf - bpb;
rnorm0 = norm(res);
disp(['|res0|=',num2str(rnorm0)]);

for cycle = 1:100
    
    % poission-boltzmann source term
    Apb = kappa2.*cosh(sol).*bfrac;
    Apb = spdiags(Apb,0, nx*ny,nx*ny);
    bpb = kappa2.*sinh(sol).*bfrac - kappa2.*cosh(sol).*bfrac.*sol;
    
    
    A = ALap - Apb;
    rhs = -rLap + bsurf + bpb;
    
    sol = A \ rhs;
    
    % check residual for nonlinear eq
    bpb = kappa2.*sinh(sol).*bfrac;
    res = ALap*sol + rLap - bsurf - bpb;
    rnorm = norm(res);
    disp(['cycle=',int2str(cycle), ', |res|=',num2str(rnorm)]);
    if rnorm <= rnorm0*1e-8
        disp(['converged']);
        break;
    end
end

sol = reshape(sol,nx,ny);


figure;
PBTestPlotSol;

PBTestCheckSum;
if (1)
    % check numerical flux through internal particles
    % they should be equal to the numerical flux through external wall
    flux_part = 0;
    for j = 2:ny-1
    for i = 2:nx-1
    % if frac(i,j)>0 | frac(i-1,j)>0 | frac(i+1,j)>0 | frac(i,j-1)>0 | frac(i,j+1)>0
    if frac(i,j)>0 | dfrac(i,j)~=0
    % if frac(i,j)>0 
    % if dfrac(i,j)>0 
        ff = 0;
        ff = ff + (sol(i+1,j)-sol(i,j)) / dx * dy;
        ff = ff + (sol(i-1,j)-sol(i,j)) / dx * dy;
        ff = ff + (sol(i,j+1)-sol(i,j)) / dy * dx;
        ff = ff + (sol(i,j-1)-sol(i,j)) / dy * dx;
        flux_part = flux_part + ff;
        % vv = kappa2*sinh(sol(i,j)) * dx*dy*(1-frac(i,j));
        % flux_part = flux_part - vv;
    end
    end
    end
    flux_part
end











