
% clc;
clear all;

% domain size and physical parameters
Lx = 1.0;
Ly = 1.0;

gx = 0.0;
gy = -100.0;

% density
rho1 = 1.0;
rho2 = 2.0;
% viscosity
m0 = 0.01;
rro = rho1;

unorth = 0;
usouth = 0;
veast = 0;
vwest = 0;

time = 0.0;

% drop size and location
rad = 0.15;
xc = 0.5;
yc = 0.7;

% numerical parameters
nx = 32;
ny = 32;
dt = 0.00125;
nstep = 100;
% nstep = 200;
% nstep = 250;
maxit = 200;
maxerr = 0.001;
beta = 1.2;

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



% set the drop
r = zeros(nx+2,ny+2) + rho1;
for i = 2:nx+1
for j = 2:ny+1
    if (x(i)-xc)^2 + (y(j)-yc)^2 < rad^2
        r(i,j) = rho2;
    end
end
end

% setup front elements
Nf = 100;
xf = zeros(1,Nf+2);
yf = zeros(1,Nf+2);
uf = zeros(1,Nf+2);
vf = zeros(1,Nf+2);
tx = zeros(1,Nf+2);
ty = zeros(1,Nf+2);

for l = 1:Nf+2
    xf(l) = xc - rad*sin(2.0*pi*(l-1)/Nf);
    yf(l) = yc + rad*cos(2.0*pi*(l-1)/Nf);
    % xf(l) = xc + rad*cos(2.0*pi*(l-1)/Nf);
    % yf(l) = yc + rad*sin(2.0*pi*(l-1)/Nf);
end

% main loop
for is = 1:nstep
    disp(['step=',int2str(is)]);
    
    %
    fx = zeros(nx+2,ny+2);
    fy = zeros(nx+2,ny+2);
    
    % fill boundary velocity
    u(1:nx+1,1) = 2*usouth - u(1:nx+1,2);
    u(1:nx+1,ny+2) = 2*unorth - u(1:nx+1,ny+1);
    v(1,1:ny+1) = 2*vwest - v(2,1:ny+1);
    v(nx+2,1:ny+1) = 2*veast - v(nx+1,1:ny+1);
    
    % U predictor
    for i = 2:nx
    for j = 2:ny+1
        ut(i,j) = u(i,j) + dt * (-0.25 * ( ...
        ((u(i+1,j)+u(i,j))^2 - (u(i,j)+u(i-1,j))^2) / dx + ...
        ((u(i,j+1)+u(i,j))*(v(i+1,j)+v(i,j)) - (u(i,j)+u(i,j-1))*(v(i+1,j-1)+v(i,j-1))) / dy) + ...
        (m0 / (0.5*(r(i+1,j)+r(i,j)))) * ( ...
        (u(i+1,j) - 2*u(i,j) + u(i-1,j)) / dx^2 + ...
        (u(i,j+1) - 2*u(i,j) + u(i,j-1)) / dy^2) + gx);
    end
    end
    
    % V predictor
    for i = 2:nx+1
    for j = 2:ny
        vt(i,j) = v(i,j) + dt * (-0.25 * ( ...
        ((u(i,j+1)+u(i,j))*(v(i+1,j)+v(i,j)) - (u(i-1,j+1)+u(i-1,j))*(v(i,j)+v(i-1,j))) / dx + ...
        ((v(i,j+1)+v(i,j))^2 - (v(i,j)+v(i,j-1))^2) / dy) + ...
        (m0 / (0.5*(r(i,j+1)+r(i,j)))) * ( ...
        (v(i+1,j) - 2*v(i,j) + v(i-1,j)) / dx^2 + ...
        (v(i,j+1) - 2*v(i,j) + v(i,j-1)) / dy^2) + gy);
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
            disp('Pressure convergence'); break;
        end
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
        ip = floor(xf(l)/dx) + 1;
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
    % move front
    for l = 2:Nf+1
        xf(l) = xf(l) + dt*uf(l);
        yf(l) = yf(l) + dt*vf(l);
    end
    xf(1) = xf(Nf+1);
    yf(1) = yf(Nf+1);
    xf(Nf+2) = xf(2);
    yf(Nf+2) = yf(2);
    
    % add points to the front
    xfold = xf;
    yfold = yf;
    j = 1;
    for l = 2:Nf+1
        ds = sqrt(((xfold(l)-xf(j))/dx)^2 + ((yfold(l)-yf(j))/dy)^2);
        if (ds > 0.5)
            j = j + 1;
            xf(j) = 0.5 * (xfold(l)+xf(j-1));
            yf(j) = 0.5 * (yfold(l)+yf(j-1));
            j = j + 1;
            xf(j) = xfold(l);
            yf(j) = yfold(l);
        elseif (ds < 0.25)
            % do nothing
        else
            j = j + 1;
            xf(j) = xfold(l);
            yf(j) = yfold(l);
        end
    end
    Nf = j - 1;
    xf(1) = xf(Nf+1);
    yf(1) = yf(Nf+1);
    xf(Nf+2) = xf(2);
    yf(Nf+2) = yf(2);
    
    % distribute gradient
    fx = zeros(nx+2,ny+2);
    fy = zeros(nx+2,ny+2);
    
    for l = 2:Nf+1
        % normal vector
        nfx = -0.5 * (yf(l+1)-yf(l-1)) * (rho2-rho1);
        nfy = 0.5 * (xf(l+1)-xf(l-1)) * (rho2-rho1);
        
        % UMAC position
        ip = floor(xf(l)/dx) + 1;
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
    
    frhs = zeros(nx+2,ny+2);
    frhs(2:nx+1,2:ny+1) = dx*fx(1:nx,2:ny+1) - dx*fx(2:nx+1,2:ny+1) + ...
        dy*fy(2:nx+1,1:ny) - dy*fy(2:nx+1,2:ny+1);
    
    if (0)
        nband = 3;
        mask = zeros(nx+2,ny+2);
        % cell position
        for l = 2:Nf+1
            ip = floor(xf(l) / dx) + 2;
            jp = floor(yf(l) / dy) + 2;
            irange = max(ip-nband,2) : min(ip+nband,nx+1);
            jrange = max(jp-nband,2) : min(jp+nband,ny+1);
            mask(irange,jrange) = 1;
        end
    end
    
    % construct density
    for iter = 1:maxit
        rsave = r;
        for i = 2:nx+1
        for j = 2:ny+1
            if (1)
            % r(i,j) = 0.25 * (r(i+1,j) + r(i-1,j) + r(i,j+1) + r(i,j-1) + ...
            % dx*fx(i-1,j) - dx*fx(i,j) + dy*fy(i,j-1) - dy*fy(i,j));
                r(i,j) = 0.25 * (r(i+1,j) + r(i-1,j) + r(i,j+1) + r(i,j-1) + frhs(i,j));
            else
                if mask(i,j) == 1
                r(i,j) = 0.25 * (r(i+1,j) + r(i-1,j) + r(i,j+1) + r(i,j-1) + frhs(i,j));
                end
            end
        end
        end
        if max(max(abs(rsave-r))) < maxerr
            disp('Density convergence'); break;
        end
    end
    
    time = time + dt;
    
    % plot results
    uu(1:nx+1,1:ny+1) = 0.5 * (u(1:nx+1,2:ny+2) + u(1:nx+1,1:ny+1));
    vv(1:nx+1,1:ny+1) = 0.5 * (v(2:nx+2,1:ny+1) + v(1:nx+1,1:ny+1));
    
    % hold off;
    % hold on;
    % imagesc(x,y,mask');
    % contourf(x,y,mask');
    
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
    drawnow
end % end main loop
























