
clear all;

%%
%% 1st-order sharp-interface like Immersed Interface method.
%%
%% For BC on immersed interface, the governing equation is extended
%% to the entire domain with augment by special jump conditions.
%% Dirichlet type (BC=1): [u]=0, [u_n]=v (unknown)
%% Neumann type (BC=0): [u]=v (unknown), [u_n]=0
%% The linear system is well-conditioned and can reuse fast Poisson solver.
%%
%% NOTE for Neumann BC, u_n is usually some 'density' in physical sense, e.g. charge.
%% In the current implementation, we may need scaling to conserve the total 'density'
%% by a factor of pi/(2*Ndim)
%%


TestHelmholtzSetup2D;

%
nsegmax = max(nx,ny)*10;
segdir = zeros(nsegmax,1);
segvec = zeros(nsegmax,2);
segidx = zeros(nsegmax,2);
segidxp = zeros(nsegmax,2);
segidxm = zeros(nsegmax,2);
seghf = zeros(nsegmax,1);

hfmin = 1.0e-3;
hfmax = 1.0-hfmin;
nseg = 0;
for j = 2:ny-1
for i = 2:nx-1
    if tag(i,j) == 0
        if tag(i-1,j) == 1
            nseg = nseg + 1;
            segdir(nseg) = 1;
            segvec(nseg,:) = [ -1,0 ];
            segidx(nseg,:) = [ i,j ];
            seghf(nseg) = sdf(i-1,j)/(sdf(i-1,j)-sdf(i,j));
            segidxp(nseg,:) = [ i-1,j ];
            segidxm(nseg,:) = [ i,j ];
        end
        if tag(i+1,j) == 1
            nseg = nseg + 1;
            segdir(nseg) = 1;
            segvec(nseg,:) = [ 1,0 ];
            segidx(nseg,:) = [ i+1,j ];
            seghf(nseg) = sdf(i+1,j)/(sdf(i+1,j)-sdf(i,j));
            segidxp(nseg,:) = [ i+1,j ];
            segidxm(nseg,:) = [ i,j ];
        end
        if tag(i,j-1) == 1
            nseg = nseg + 1;
            segdir(nseg) = 2;
            segvec(nseg,:) = [ 0,-1 ];
            segidx(nseg,:) = [ i,j ];
            seghf(nseg) = sdf(i,j-1)/(sdf(i,j-1)-sdf(i,j));
            segidxp(nseg,:) = [ i,j-1 ];
            segidxm(nseg,:) = [ i,j ];
        end
        if tag(i,j+1) == 1
            nseg = nseg + 1;
            segdir(nseg) = 2;
            segvec(nseg,:) = [ 0,1 ];
            segidx(nseg,:) = [ i,j+1 ];
            seghf(nseg) = sdf(i,j+1)/(sdf(i,j+1)-sdf(i,j));
            segidxp(nseg,:) = [ i,j+1 ];
            segidxm(nseg,:) = [ i,j ];
        end
    end
end
end
% pack segment data
segdir = segdir(1:nseg);
segvec = segvec(1:nseg,:);
segidx = segidx(1:nseg,:);
segidxp = segidxp(1:nseg,:);
segidxm = segidxm(1:nseg,:);
seghf = seghf(1:nseg);
% clip height fraction
seghf = min(seghf,hfmax);
seghf = max(seghf,hfmin);

% seghf(:) = 0.5;

disp(['Build A-Laplacian']);
tic;
[ ALap, rLap ] = MakeLap2Da(nx,ny,dx,dy, bctype,bcval);
ALap = ALap - kappa2.*speye(ncell);
b1 = -rLap;
% asolve = @(b) SolveIChol(b,RLap,RLapt,[]);
asolve = @(b) ALap \ b;
toc;

disp(['Build B-matrix']);
tic;
Bmat = sparse(ncell,nseg);
for iseg = 1:nseg
    i = segidx(iseg,1);
    j = segidx(iseg,2);
    dh = celldh(segdir(iseg));
    ip = segidxp(iseg,1);
    jp = segidxp(iseg,2);
    im = segidxm(iseg,1);
    jm = segidxm(iseg,2);
    indp = ip + (jp-1)*nx;
    indm = im + (jm-1)*nx;
    theta = seghf(iseg);
    
    ubc = partbc(owner(im,jm));
    if ubc == 0
        Bmat(indp,iseg) = 1.0/dh^2;
        Bmat(indm,iseg) = -1.0/dh^2;
    elseif ubc == 1
        Bmat(indp,iseg) = -(1-theta)/dh;
        Bmat(indm,iseg) = -theta/dh;
    end
end
toc;

disp('Build C-matrix');
tic;
Cmat = sparse(nseg,ncell);
for iseg = 1:nseg
    i = segidx(iseg,1);
    j = segidx(iseg,2);
    dh = celldh(segdir(iseg));
    ip = segidxp(iseg,1);
    jp = segidxp(iseg,2);
    im = segidxm(iseg,1);
    jm = segidxm(iseg,2);
    
    theta = seghf(iseg);
    indp = ip + (jp-1)*nx;
    indm = im + (jm-1)*nx;
    
    ubc = partbc(owner(im,jm));
    if ubc == 0
        Cmat(iseg,indp) = 1.0/dh^2;
        Cmat(iseg,indm) = -1.0/dh^2;
    elseif ubc == 1
        Cmat(iseg,indp) = (1-theta) / dh;
        Cmat(iseg,indm) = theta / dh;
    end
end
toc;

disp('Build D-matrix');
tic;
Dmat = sparse(nseg,nseg);
b2 = zeros(nseg,1);
for iseg = 1:nseg
    i = segidx(iseg,1);
    j = segidx(iseg,2);
    dh = celldh(segdir(iseg));
    im = segidxm(iseg,1);
    jm = segidxm(iseg,2);
    
    theta = seghf(iseg);
    ubc = partbc(owner(im,jm));
    up = partu(owner(im,jm));
    if ubc == 0
        up = up * pi/4;
        Dmat(iseg,iseg) = -1.0/dh^2;
        b2(iseg) = up/dh;
    elseif ubc == 1
        Dmat(iseg,iseg) = -theta*(1-theta);
        b2(iseg) = up/dh;
    end
end
toc;

disp('Solver');
if (0)
    mat = [ ALap, Bmat; Cmat, Dmat ];
    rhs = [ b1; b2];
    sol = mat \ rhs;
    % [sol,ret,res,iter] = bicgstab(mat,rhs,1.0e-9,1000);
    % disp(['Solver: ret=',int2str(ret),';res=',num2str(res),';iter=',int2str(iter)]);
    u = sol(1:ncell);
    u = reshape(u,nx,ny);
end
if (1)
    % afun = @(v) Dmat*v - Cmat * (ALap \ (Bmat*v));
    % vrhs = b2 - Cmat*(ALap \ b1);
    afun = @(v) Dmat*v - Cmat*asolve(Bmat*v);
    vrhs = b2 - Cmat*asolve(b1);
    tol = 1.0e-8;
    maxiter = 200;
    [vsol,ret,res,iter] = bicgstab(afun,vrhs,tol,maxiter);
    % [vsol,ret,res,iter] = gmres(afun,vrhs,10,tol,maxiter);
    disp(['Solver: ret=',int2str(ret),';res=',num2str(res),';iter=',int2str(iter)]);
    
    urhs = b1 - Bmat*vsol;
    usol = ALap \ urhs;
    u = reshape(usol,nx,ny);
end
if (0)
    umat = ALap - Bmat*inv(Dmat)*Cmat;
    urhs = -rLap - Bmat*inv(Dmat)*b2;
    usol = umat \ urhs;
    u = reshape(usol,nx,ny);
end

u(tag_nonfd) = nan;
if (1)
    figure;
    hold on;
    
    % contourf(xcell,ycell,u);
    usc = imagesc([xlo+dx/2,xhi-dx/2],[ylo+dy/2,yhi-dy/2],u');
    set(usc,'AlphaData',~isnan(u'));
    % pcolor(xcell,ycell,u);
    set(gca,'YDir','normal');
    colorbar;
    
    %% draw interface shape
    contour(xcell,ycell,sdf,[0,0]);
    % mesh(nodex,nodey,zeros(nx+1,ny+1));
    
    % draw internal faces
    for iseg = 1:nseg
        ii = segidx(iseg,1);
        jj = segidx(iseg,2);
        x1 = (ii-1)*dx + xlo;
        y1 = (jj-1)*dy + ylo;
        x2 = x1; if segdir(iseg)==2; x2 = x2 + dx; end;
        y2 = y1; if segdir(iseg)==1; y2 = y2 + dy; end;
        plot([x1,x2],[y1,y2],'.-k');
    end
    
    hold off;
    axis equal;
    axis([xlo xhi ylo yhi]);
end

if (1)
    figure;
    sol = u;
    PBTestPlotSol;
end


