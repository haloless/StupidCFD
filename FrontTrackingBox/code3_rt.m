
% clc;
clear all;

% domain size 
Lx = 2.0 * pi;
Ly = 6.0 * pi;
% gravity
gx = 0.0;
gy = 1.0;
% density
rho1 = 0.25;
rho2 = 1.0;
% viscosity
mu1 = 0.00;
mu2 = 0.00;
% surface tension
sigma = 0.075;
%
rro = 0.0;

% BC velocity
% unorth = 0;
% usouth = 0;
% veast = 0;
% vwest = 0;

time = 0.0;

% RT perturbation
% delta = 0.25;
delta = 0.1 * Lx;

% numerical parameters
nx = 40;
ny = 120;
dt = 0.01;
nstep = 1000;
maxit = 1000;
maxerr = 0.001;
beta = 1.25;
% beta = 1.5;

% allocate buffers
u = zeros(nx+1,ny+2);
v = zeros(nx+2,ny+1);
p = zeros(nx+2,ny+2);
ut = zeros(nx+1,ny+2);
vt = zeros(nx+2,ny+1);
uu = zeros(nx+1,ny+1);
vv = zeros(nx+1,ny+1);
tmp1 = zeros(nx+2,ny+2);
tmp2 = zeros(nx+2,ny+2);

fx = zeros(nx+2,ny+2);
fy = zeros(nx+2,ny+2);

% setup grid 
dx = Lx / nx;
dy = Ly / ny;
% cell position
for i = 1:nx+2
    x(i) = dx * (i-1.5);
end
for j = 1:ny+2
    y(j) = dy * (j-1.5);
end
% edge position
for i = 1:nx+1
    xh(i) = dx*(i-1);
end
for j = 1:ny+1
    yh(j) = dy*(j-1);
end



% set density and viscosity
r = zeros(nx+2,ny+2);
m = zeros(nx+2,ny+2);
for i = 2:nx+1
for j = 2:ny+1
    if y(j) >= Ly*0.5+delta*cos(x(i))
        r(i,j) = rho2;
        m(i,j) = mu2;
    else
        r(i,j) = rho1;
        m(i,j) = mu1;
    end
end
end
r(1,:) = r(2,:); r(nx+2,:) = r(nx+1,:); r(:,1) = r(:,2); r(:,ny+2) = r(:,ny+1);
m(1,:) = m(2,:); m(nx+2,:) = m(nx+1,:); m(:,1) = m(:,2); m(:,ny+2) = m(:,ny+1);

% setup front elements
Nf = 100;
xf = zeros(1,Nf+2);
yf = zeros(1,Nf+2);
uf = zeros(1,Nf+2);
vf = zeros(1,Nf+2);
tx = zeros(1,Nf+2);
ty = zeros(1,Nf+2);

for l = 1:Nf+2
    xf(l) = Lx/(Nf-1) * (l-2);
    yf(l) = Ly/2 + delta*cos(xf(l));
end

if (0)
    uu(1:nx+1,1:ny+1) = 0.5 * (u(1:nx+1,2:ny+2) + u(1:nx+1,1:ny+1));
    vv(1:nx+1,1:ny+1) = 0.5 * (v(2:nx+2,1:ny+1) + v(1:nx+1,1:ny+1));
    
    % hold off;
    % hold on;
    % imagesc(x,y,mask');
    % contourf(x,y,mask');
    
    contour(x,y,flipud(rot90(r)));
    % contour(x,y,r',[0.5*(rho2+rho1),0.5*(rho2+rho1)],'linewidth',3);
    % axis equal;
    axis([0 Lx 0 Ly]);
    hold on;
    quiver(xh,yh,flipud(rot90(uu)),flipud(rot90(vv)),'r');
    plot(xf(2:Nf+1),yf(2:Nf+1),'k.-','linewidth',1); 
    % pause(0.01);
    
    % title(['step=',int2str(is),';time=',num2str(time)]);
    hold off;
    drawnow
    
    return
end

% main loop
for is = 1:nstep
    disp(['step=',int2str(is)]);
        
    % surface tension
    % find tangent vectors
    for l = 1:Nf+1
        dxf = xf(l+1) - xf(l);
        dyf = yf(l+1) - yf(l);
        ds = sqrt(dxf^2 + dyf^2);
        tx(l) = dxf / ds;
        ty(l) = dyf / ds;
    end
    % tx(Nf+2) = tx(2);
    % ty(Nf+2) = ty(2);
    
    % distribute to the grid
    fx = zeros(nx+2,ny+2);
    fy = zeros(nx+2,ny+2);
    for l = 2:Nf+1
        nfx = sigma * (tx(l)-tx(l-1));
        nfy = sigma * (ty(l)-ty(l-1));
        
        % UMAC position
        ip = floor(xf(l)/dx) + 1; ip = min(ip,nx); ip = max(ip,1);
        jp = floor((yf(l)+0.5*dy)/dy) + 1;
        ax = xf(l)/dx - ip + 1;
        ay = (yf(l)+0.5*dy)/dy - jp + 1;
        fx(ip,jp) = fx(ip,jp) + (1.0-ax)*(1.0-ay)*nfx/dx/dy;
        fx(ip+1,jp) = fx(ip+1,jp) + ax*(1.0-ay)*nfx/dx/dy;
        fx(ip,jp+1) = fx(ip,jp+1) + (1.0-ax)*ay*nfx/dx/dy;
        fx(ip+1,jp+1) = fx(ip+1,jp+1) + ax*ay*nfx/dx/dy;
        
        % VMAC position
        ip = floor((xf(l)+0.5*dx)/dx) + 1;
        jp = floor(yf(l)/dy) + 1;
        ax = (xf(l)+0.5*dx)/dx - ip + 1;
        ay = yf(l)/dy - jp + 1;
        fy(ip,jp) = fy(ip,jp) + (1.0-ax)*(1.0-ay)*nfy/dx/dy;
        fy(ip+1,jp) = fy(ip+1,jp) + ax*(1.0-ay)*nfy/dx/dy;
        fy(ip,jp+1) = fy(ip,jp+1) + (1.0-ax)*ay*nfy/dx/dy;
        fy(ip+1,jp+1) = fy(ip+1,jp+1) + ax*ay*nfy/dx/dy;
    end
    
    
    % CFD
    % fill boundary velocity
    % u(1:nx+1,1) = 2*usouth - u(1:nx+1,2);
    % u(1:nx+1,ny+2) = 2*unorth - u(1:nx+1,ny+1);
    % v(1,1:ny+1) = 2*vwest - v(2,1:ny+1);
    % v(nx+2,1:ny+1) = 2*veast - v(nx+1,1:ny+1);
    % free-slip BC
    u(1:nx+1,1) = u(1:nx+1,2);
    u(1:nx+1,ny+2) = u(1:nx+1,ny+1);
    v(1,1:ny+1) = v(2,1:ny+1);
    v(nx+2,1:ny+1) = v(nx+1,1:ny+1);
    
    % U advection
    for i = 2:nx
    for j = 2:ny+1
        ue = 0.5 * (u(i+1,j) + u(i,j));
        uw = 0.5 * (u(i,j) + u(i-1,j));
        unn = 0.5 * (u(i,j+1) + u(i,j));
        us = 0.5 * (u(i,j) + u(i,j-1));
        vnn = 0.5 * (v(i,j) + v(i+1,j));
        vs = 0.5 * (v(i,j-1) + v(i+1,j-1));
        
        rp = 0.5 * (r(i+1,j) + r(i,j));
        
        ut(i,j) = u(i,j) + dt * ( ...
            -((ue^2-uw^2)/dx + (unn*vnn-us*vs)/dy) ... 
            + fx(i,j) / rp ...
            - (1.0 - rro/rp) * gx);
    end
    end
    
    % V advection
    for i = 2:nx+1
    for j = 2:ny
        ve = 0.5 * (v(i+1,j) + v(i,j));
        vw = 0.5 * (v(i,j) + v(i-1,j));
        vnn = 0.5 * (v(i,j+1) + v(i,j));
        vs = 0.5 * (v(i,j) + v(i,j-1));
        ue = 0.5 * (u(i,j) + u(i,j+1));
        uw = 0.5 * (u(i-1,j) + u(i-1,j+1));
        
        rp = 0.5 * (r(i,j+1) + r(i,j));
        
        vt(i,j) = v(i,j) + dt * ( ...
            -((ue*ve-uw*vw)/dx + (vnn^2-vs^2)/dy) ...
            + fy(i,j) / rp ...
            - (1.0 - rro/rp) * gy);
    end
    end
    
    % U diffusion
    for i = 2:nx
    for j = 2:ny+1
        due = m(i+1,j) * (u(i+1,j)-u(i,j)) / dx;
        duw = m(i,j) * (u(i,j)-u(i-1,j)) / dx;
        mnorth = 0.25 * (m(i,j)+m(i+1,j)+m(i+1,j+1)+m(i,j+1));
        dun = mnorth * ((u(i,j+1)-u(i,j))/dy + (v(i+1,j)-v(i,j))/dx);
        ms = 0.25 * (m(i,j)+m(i+1,j)+m(i+1,j-1)+m(i,j-1));
        dus = ms * ((u(i,j)-u(i,j-1))/dy + (v(i+1,j-1)-v(i,j-1))/dx);
        
        rp = 0.5 * (r(i+1,j) + r(i,j));
        
        ut(i,j) = ut(i,j) + dt * ( ...
            1.0/dx * 2.0 * (due - duw) + ...
            1.0/dy * (dun - dus)) / rp;
    end
    end
    
    % V diffusion
    for i = 2:nx+1
    for j = 2:ny
        me = 0.25 * (m(i,j) + m(i+1,j) + m(i+1,j+1) + m(i,j+1));
        dve = me * ((u(i,j+1)-u(i,j))/dy + (v(i+1,j)-v(i,j))/dx);
        mw = 0.25 * (m(i,j) + m(i,j+1) + m(i-1,j+1) + m(i-1,j));
        dvw = mw * ((u(i-1,j+1)-u(i-1,j))/dy + (v(i,j)-v(i-1,j))/dx);
        dvn = m(i,j+1) * (v(i,j+1)-v(i,j)) / dy;
        dvs = m(i,j) * (v(i,j)-v(i,j-1)) / dy;
        
        rp = 0.5 * (r(i,j+1) + r(i,j));
        
        vt(i,j) = vt(i,j) + dt * ( ...
            1.0/dx * (dve - dvw) + ...
            1.0/dy * 2.0 * (dvn - dvs)) / rp;
    end
    end
    
    % PPE matrix and RHS
    rt = r;
    lrg = 1000;
    rt(1:nx+2,1) = lrg;
    rt(1:nx+2,ny+2) = lrg;
    rt(1,1:ny+2) = lrg;
    rt(nx+2,1:ny+2) = lrg;
    
    for i = 2:nx+1
    for j = 2:ny+1
        tmp1(i,j) = 0.5/dt * ((ut(i,j)-ut(i-1,j))/dx + (vt(i,j)-vt(i,j-1))/dy);
        tmp2(i,j) = 1.0 / ( ...
        (1.0/dx) * (1.0/(dx*(rt(i+1,j)+rt(i,j))) + 1.0/(dx*(rt(i-1,j)+rt(i,j)))) + ...
        (1.0/dy) * (1.0/(dy*(rt(i,j+1)+rt(i,j))) + 1.0/(dy*(rt(i,j-1)+rt(i,j)))));
    end
    end
    
    % solve PPE
    if (0)
    for it = 1:maxit
        psave = p;
        for i = 2:nx+1
        for j = 2:ny+1
            p(i,j) = (1.0-beta)*p(i,j) + beta*tmp2(i,j) * ( ...
            (1.0/dx) * (p(i+1,j)/(dx*(rt(i+1,j)+rt(i,j))) + p(i-1,j)/(dx*(rt(i-1,j)+rt(i,j)))) + ...
            (1.0/dy) * (p(i,j+1)/(dy*(rt(i,j+1)+rt(i,j))) + p(i,j-1)/(dy*(rt(i,j-1)+rt(i,j)))) - ...
            tmp1(i,j));
        end
        end
        if max(max(abs(psave-p))) < maxerr
            disp(['Pressure convergence, iter=',int2str(it)]); break;
        end
    end
    else
        Lap = sparse(nx*ny,nx*ny);
        idx = 1;
        for j = 2:ny+1
        for i = 2:nx+1
            if j ~= 2
                coef = 1.0 / dy^2 / (rt(i,j)+rt(i,j-1));
                Lap(idx,idx-nx) = coef;
                Lap(idx,idx) = Lap(idx,idx) - coef;
            end
            if i ~= 2
                coef = 1.0 / dx^2 / (rt(i,j)+rt(i-1,j));
                Lap(idx,idx-1) = coef;
                Lap(idx,idx) = Lap(idx,idx) - coef;
            end
            if i ~= nx+1
                coef = 1.0 / dx^2 / (rt(i,j)+rt(i+1,j));
                Lap(idx,idx+1) = coef;
                Lap(idx,idx) = Lap(idx,idx) - coef;
            end
            if j ~= ny+1
                coef = 1.0/dy^2 / (rt(i,j)+rt(i,j+1));
                Lap(idx,idx+nx) = coef;
                Lap(idx,idx) = Lap(idx,idx) - coef;
            end
            idx = idx + 1;
        end
        end
        rhs = reshape(tmp1(2:nx+1,2:ny+1), nx*ny,1);
        Lap(1,:) = 0; Lap(:,1) = 0; Lap(1,1) = 1;
        rhs(1) = 0;
        
        sol = Lap \ rhs;
        p(2:nx+1,2:ny+1) = reshape(sol, nx,ny);
        
    end
    
    % U corrector
    for i = 2:nx
    for j = 2:ny+1
        u(i,j) = ut(i,j) - dt*2.0/dx * (p(i+1,j)-p(i,j)) / (r(i+1,j)+r(i,j));
    end
    end
    % V corrector
    for i = 2:nx+1
    for j = 2:ny
        v(i,j) = vt(i,j) - dt*2.0/dy * (p(i,j+1)-p(i,j)) / (r(i,j+1)+r(i,j));
    end
    end
    
    % advect front
    % front velocity
    for l = 2:Nf+1
        ip = floor(xf(l)/dx) + 1; ip = min(ip,nx); ip = max(ip,1);
        jp = floor((yf(l)+0.5*dy)/dy) + 1;
        ax = (xf(l))/dx - ip + 1;
        ay = (yf(l)+0.5*dy)/dy - jp + 1;
        uf(l) = (1.0-ax)*(1.0-ay)*u(ip,jp) + ax*(1.0-ay)*u(ip+1,jp) + ...
        (1.0-ax)*ay*u(ip,jp+1) + ax*ay*u(ip+1,jp+1);
        
        ip = floor((xf(l)+0.5*dx)/dx) + 1;
        jp = floor(yf(l)/dy) + 1;
        ax = (xf(l)+0.5*dx)/dx - ip + 1;
        ay = (yf(l))/dy - jp + 1;
        vf(l) = (1.0-ax)*(1.0-ay)*v(ip,jp) + ax*(1.0-ay)*v(ip+1,jp) + ...
        (1.0-ax)*ay*v(ip,jp+1) + ax*ay*v(ip+1,jp+1);
    end
    % constraint for nodes attached to wall
    uf(2) = 0;
    uf(Nf+1) = 0;
    % move front
    for l = 2:Nf+1
        xf(l) = xf(l) + dt*uf(l);
        yf(l) = yf(l) + dt*vf(l);
    end
    % xf(2) = 0;
    % xf(Nf+1) = Lx;
    xf(1) = -xf(3);
    yf(1) = yf(3);
    xf(Nf+2) = Lx*2.0 - xf(Nf);
    yf(Nf+2) = yf(Nf);
    
    % add points to the front
    xfold = xf;
    yfold = yf;
    j = 2;
    for l = 3:Nf+1
        ds = sqrt(((xfold(l)-xf(j))/dx)^2 + ((yfold(l)-yf(j))/dy)^2);
        if (ds > 0.5)
            j = j + 1;
            xf(j) = 0.5 * (xfold(l)+xf(j-1));
            yf(j) = 0.5 * (yfold(l)+yf(j-1));
            j = j + 1;
            xf(j) = xfold(l);
            yf(j) = yfold(l);
        % elseif (ds < 0.25)
            % do nothing
        else
            j = j + 1;
            xf(j) = xfold(l);
            yf(j) = yfold(l);
        end
    end
    Nf = j - 1;
    % xf(1) = xf(Nf+1);
    % yf(1) = yf(Nf+1);
    % xf(Nf+2) = xf(2);
    % yf(Nf+2) = yf(2);
    % xf(2) = 0;
    % xf(Nf+1) = Lx;
    xf(1) = -xf(3);
    yf(1) = yf(3);
    xf(Nf+2) = Lx*2.0 - xf(Nf);
    yf(Nf+2) = yf(Nf);

    
    % distribute gradient
    fx = zeros(nx+2,ny+2);
    fy = zeros(nx+2,ny+2);
    
    for l = 2:Nf+1
        % normal vector
        nfx = -0.5 * (yf(l+1)-yf(l-1)) * (rho2-rho1);
        nfy = 0.5 * (xf(l+1)-xf(l-1)) * (rho2-rho1);
        
        % UMAC position
        ip = floor(xf(l)/dx) + 1; ip = min(ip,nx); ip = max(ip,1);
        jp = floor((yf(l)+0.5*dy)/dy) + 1;
        ax = xf(l)/dx - ip + 1;
        ay = (yf(l)+0.5*dy)/dy - jp + 1;
        fx(ip,jp) = fx(ip,jp) + (1.0-ax)*(1.0-ay)*nfx/dx/dy;
        fx(ip+1,jp) = fx(ip+1,jp) + ax*(1.0-ay)*nfx/dx/dy;
        fx(ip,jp+1) = fx(ip,jp+1) + (1.0-ax)*ay*nfx/dx/dy;
        fx(ip+1,jp+1) = fx(ip+1,jp+1) + ax*ay*nfx/dx/dy;
        
        % VMAC position
        ip = floor((xf(l)+0.5*dx)/dx) + 1;
        jp = floor(yf(l)/dy) + 1;
        ax = (xf(l)+0.5*dx)/dx - ip + 1;
        ay = yf(l)/dy - jp + 1;
        fy(ip,jp) = fy(ip,jp) + (1.0-ax)*(1.0-ay)*nfy/dx/dy;
        fy(ip+1,jp) = fy(ip+1,jp) + ax*(1.0-ay)*nfy/dx/dy;
        fy(ip,jp+1) = fy(ip,jp+1) + (1.0-ax)*ay*nfy/dx/dy;
        fy(ip+1,jp+1) = fy(ip+1,jp+1) + ax*ay*nfy/dx/dy;
    end
    
    rrhs = zeros(nx+2,ny+2);
    rrhs(2:nx+1,2:ny+1) = dx*fx(1:nx,2:ny+1) - dx*fx(2:nx+1,2:ny+1) + ...
        dy*fy(2:nx+1,1:ny) - dy*fy(2:nx+1,2:ny+1);
    
    % construct density
    for iter = 1:maxit
        rsave = r;
        for i = 2:nx+1
        for j = 2:ny+1
            r(i,j) = (1.0-beta)*r(i,j) + beta * ...
                0.25 * (r(i+1,j) + r(i-1,j) + r(i,j+1) + r(i,j-1) + rrhs(i,j));
        end
        end
        r(1,:) = r(2,:); r(nx+2,:) = r(nx+1,:); r(:,1) = r(:,2); r(:,ny+2) = r(:,ny+1);
        if max(max(abs(rsave-r))) < maxerr
            disp(['Density convergence, iter=',int2str(iter)]); break;
        end
    end
    
    % update viscosity
    m = zeros(nx+2,ny+2) + mu1;
    for i = 2:nx+1
    for j = 2:ny+1
        m(i,j) = mu1 + (mu2-mu1) * (r(i,j)-rho1)/(rho2-rho1);
    end
    end
    m(1,:) = m(2,:); m(nx+2,:) = m(nx+1,:); m(:,1) = m(:,2); m(:,ny+2) = m(:,ny+1);
    
    
    % if (0)
        % nband = 3;
        % mask = zeros(nx+2,ny+2);
        % % cell position
        % for l = 2:Nf+1
            % ip = floor(xf(l) / dx) + 2;
            % jp = floor(yf(l) / dy) + 2;
            % irange = max(ip-nband,2) : min(ip+nband,nx+1);
            % jrange = max(jp-nband,2) : min(jp+nband,ny+1);
            % mask(irange,jrange) = 1;
        % end
    % end
    

    time = time + dt;
    
    % plot results
    uu(1:nx+1,1:ny+1) = 0.5 * (u(1:nx+1,2:ny+2) + u(1:nx+1,1:ny+1));
    vv(1:nx+1,1:ny+1) = 0.5 * (v(2:nx+2,1:ny+1) + v(1:nx+1,1:ny+1));
    
    subplot(1,2,1);
    contour(x,y,flipud(rot90(r)));
    % contour(x,y,r',[0.5*(rho2+rho1),0.5*(rho2+rho1)],'linewidth',3);
    axis equal;
    axis([0 Lx 0 Ly]);
    hold on;
    quiver(xh,yh,flipud(rot90(uu)),flipud(rot90(vv)),'r');
    plot(xf(1:Nf),yf(1:Nf),'k.-','linewidth',1); 
    % pause(0.01);
    
    title(['step=',int2str(is),';time=',num2str(time)]);
    hold off;
    
    
    subplot(1,2,2);
    contourf(x(2:nx+1),y(2:ny+1),p(2:nx+1,2:ny+1)');
    title(['pressure']);
    
    hold off;
    drawnow
end % end main loop
























