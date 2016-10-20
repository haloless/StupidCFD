
clear all;

%
PNPTestSetup2D;

%
% External electric field
%
disp(['Build External Laplacian']);
[ Aext, rext ] = MakeLap2Da(nx,ny,dx,dy, bctype_ext,bcval_ext);

disp(['Modify External Laplacian']);
[ Aext, rext ] = ModifyLap2D(Aext,rext, nx,ny,dx,dy,xlo,ylo, sdf,tag,owner, npart,partbc_ext,partval_ext);

disp(['Solve External']);
rhs = -rext;
sol = zeros(nx*ny,1);
[Lext,Uext] = ilu(Aext,struct('type','nofill', 'droptol',1.0e-6));
[sol,flag,relres,iter] = bicgstab(Aext, rhs, 1e-8, 2000, Lext,Uext,sol);
disp(['solver: flag=',int2str(flag), '; res=',num2str(relres), '; iter=',int2str(iter)]);

% external potential
psi_ext = reshape(sol,nx,ny);
% its gradient
[gphix_ext,gphiy_ext] = CalcGradPhi(psi_ext, nx,ny,dx,dy, bctype_ext,bcval_ext);

% plot External potential
if (1)
    figure;
    solrange = [-0.9:0.2:-0.1,0.0,0.1:0.2:0.9];
    PNPTestPlotSol;
    title('\psi Ext');
    hold on;
    % quiver(xumac,yumac,-gphix_ext,zeros(nx+1,ny));
    % quiver(xvmac,yvmac,zeros(nx,ny+1),-gphiy_ext);
    quiver(xcell,ycell,-0.5*(gphix_ext(1:nx,:)+gphix_ext(2:nx+1,:)),-0.5*(gphiy_ext(:,1:ny)+gphiy_ext(:,2:ny+1)));
    hold off;
end

%
% Coupled internal field & charge density
%

% basic Laplacian
disp(['Build Internal Laplacian']);
[ Aint,rint ] = MakeLap2Da(nx,ny,dx,dy, bctype_int,bcval_int);

% internal forcing
disp(['Modify Internal Laplacian']);
[ Aint,rint ] = ModifyLap2D(Aint,rint, nx,ny,dx,dy,xlo,ylo, sdf,tag,owner, npart,partbc_int,partval_int);

% create preconditioner
[Lint,Uint] = ilu(Aint,struct('type','nofill', 'droptol',1.0e-6));

% relaxation
wpsi = 0.1;
wspec = 0.1;

psi_int = zeros(nx,ny);

figdbg = figure;

for cycle = 1:50
    disp(['cycle=',int2str(cycle)]);
    
    % derive charge density
    rho = zeros(nx,ny);
    for ispec = 1:nspec
        rho = rho + specz(ispec).*specdens(:,:,ispec);
    end
    rho = reshape(rho,nx*ny,1);
    
    % solve internal potential
    rhs = -clambda2 .* rho;
    rhs = rhs - rint;
    sol = reshape(psi_int,nx*ny,1);
    disp(['Solve Internal']);
    [sol,flag,relres,iter] = bicgstab(Aint, rhs, 1e-8, 2000, Lint,Uint,sol);
    disp(['solver: flag=',int2str(flag), '; res=',num2str(relres), '; iter=',int2str(iter)]);
    
    % relaxation
    sol = reshape(sol,nx,ny);
    psi_int = wpsi.*sol + (1-wpsi).*psi_int;
    
    % calculate gradient
    [gphix_int,gphiy_int] = CalcGradPhi(psi_int,nx,ny,dx,dy, bctype_int,bcval_int);
    
    gphix = gphix_ext + gphix_int;
    gphiy = gphiy_ext + gphiy_int;
    
    % solve
    for ispec = 1:nspec
        [sol,flag,relres,iter] = PNPSolveSpec(specdens(:,:,ispec),specz(ispec), ...
        gphix,gphiy, sdf,tag, nx,ny,dx,dy,bctype_int);
        disp(['solver: spec=',int2str(ispec), ...
        '; flag=',int2str(flag), '; res=',num2str(relres), '; iter=',int2str(iter)]);
        
        cnew = reshape(sol,nx,ny);
        specdens(:,:,ispec) = wspec.*cnew + (1-wspec).*specdens(:,:,ispec);
    end
    
    figure(figdbg);
    if (1)
    subplot(3,1,1);
    sol = psi_int;
    solrange = -1:0.05:0;
    PNPTestPlotSol;
    title('\psi Int');
    subplot(3,1,2);
    sol = specdens(:,:,1); 
    sol(tag_nonfd) = nan;
    solrange = [];
    PNPTestPlotSol;
    title('cation');
    subplot(3,1,3);
    sol = specdens(:,:,2);
    sol(tag_nonfd) = nan;
    solrange = [];
    PNPTestPlotSol;
    title('anion');
    end
    pause;
end

figure;
sol = psi_int;
solrange = -1:0.05:0;
PNPTestPlotSol;
title('\psi Int');
% hold on;
% quiver(xcell,ycell,-0.5*(gphix_int(1:nx,:)+gphix_int(2:nx+1,:)),-0.5*(gphiy_int(:,1:ny)+gphiy_int(:,2:ny+1)));
% hold off;

figure;
sol = specdens(:,:,1); 
sol(tag_nonfd) = nan;
solrange = [];
PNPTestPlotSol;
title('cation');

figure;
sol = specdens(:,:,2);
sol(tag_nonfd) = nan;
solrange = [];
PNPTestPlotSol;
title('anion');










