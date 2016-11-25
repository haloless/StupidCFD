
clear all;

global Emat Hmat;

twopi = pi * 2;
ImUnit = 1i;

rho = 1.0;
nu = 1.0;

Lx = pi * 2;
Ly = pi * 2;
xlo = 0;
ylo = 0;
xhi = Lx;
yhi = Ly;
xcen = (xlo+xhi) / 2;
ycen = (ylo+yhi) / 2;

% nx = 32;
% ny = 32;
nx = 64;
ny = 64;
% nx = 128;
% ny = 128;
dx = Lx / nx;
dy = Ly / ny;

ng = nx * ny;
matind = reshape(1:ng, nx,ny);

xs = linspace(0,Lx-dx,nx);
ys = linspace(0,Ly-dy,ny);
[xs,ys] = ndgrid(xs,ys);

xper = linspace(0,Lx,nx+1);
yper = linspace(0,Ly,ny+1);
[xper,yper] = ndgrid(xper,yper);

kx = zeros(nx,ny);
ky = zeros(nx,ny);
for j = 1:ny
    kx(:,j) = [0:nx/2-1, -nx/2:-1];
end
for i = 1:nx
    ky(i,:) = [0:ny/2-1, -ny/2:-1];
end
k2 = kx.^2 + ky.^2;

fs = zeros(nx,ny);
gs = zeros(nx,ny);
% fs(:,:) = sin(xs) .* cos(xs) + 2*nu*cos(xs).*sin(ys);
% gs(:,:) = sin(ys) .* cos(ys) - 2*nu*sin(xs).*cos(ys);


% hp = dx * 1.5;
% hp = dx * 1.7;
hp = dx * 1.8;
% hp = dx * 2.0;
% hp = dx * 2.1;
% hp = dx * 2.5;
% hp = dx * 3.0;

%
if 0
	%
	rad = 1.0;
	np = round(2*pi*rad / hp);
	phi = (0:np-1)' .* (2*pi/np);
	xp = rad.*cos(phi) + xcen;
	yp = rad.*sin(phi) + ycen;
	%
	up = zeros(np,1);
	vp = zeros(np,1);
	up(:) = 1.0;
end
if 0
	rout = 2.0;
	rin = 1.0;
	wout = 0;
	win = 1.0;
	
	npout = round(2*pi*rout / hp);
	npin = round(2*pi*rin / hp);
	np = npout + npin;
	
	xp = zeros(np,1);
	yp = zeros(np,1);
	up = zeros(np,1);
	vp = zeros(np,1);
	
	tmp = (0:npout-1)' .* (2*pi/npout);
	xout = rout.*cos(tmp);
	yout = rout.*sin(tmp);
	xp(1:npout) = xout + xcen;
	yp(1:npout) = yout + ycen;
	up(1:npout) = -wout .* yout;
	vp(1:npout) = wout .* xout;
	
	tmp = (0:npin-1)' .* (2*pi/npin);
	xin = rin.*cos(tmp);
	yin = rin.*sin(tmp);
	xp(npout+1:np) = xin + xcen;
	yp(npout+1:np) = yin + ycen;
	up(npout+1:np) = -win .* yin;
	vp(npout+1:np) = win .* xin;
end



if 1
	figure;
	plot(xp,yp,'.k');
	hold on;
	quiver(xp,yp,up,vp);
	hold off;
	axis equal;
	axis([0 Lx 0 Ly]);
end

% dfunc = @(r) 0.25 .* (1 + cos(0.5*pi.*r));
dfunc = @(r) 0.25 .* (1 + cos(pi.*r));
dspan = 3;
dsupp = 2;

Emat = sparse(ng,np);
for ip = 1:np
	x0 = xp(ip);
	y0 = yp(ip);
	% locate in grid
	i0 = round(x0 / dx) + 1;
	j0 = round(y0 / dy) + 1;
	%
	for ii = i0-dspan:i0+dspan
	for jj = j0-dspan:j0+dspan
		xx = abs(dx*(ii-1) - x0) / (dsupp*dx);
		yy = abs(dy*(jj-1) - y0) / (dsupp*dy);
		
		rx = 0;
		if xx < 1
			% rx = 1.0/dx * dfunc(xx);
			rx = dfunc(xx);
		end
		ry = 0;
		if yy < 1
			% ry = 1.0/dy * dfunc(yy);
			ry = dfunc(yy);
		end
		
		dd = rx * ry;
		if dd > 0
			% find real grid index
			ig = ii;
			jg = jj;
			if ig < 1
				ig = ig + nx;
			elseif ig > nx
				ig = ig - nx;
			end
			if jg < 1
				jg = jg + ny;
			elseif jg > ny
				jg = jg - ny;
			end
			
			ind = matind(ig,jg);
			Emat(ind,ip) = dd;
		end
	end
	end
end
Hmat = Emat';




rhs = [up; vp];
afun = @(pp) IBSolve2d(kx,ky, np,pp(1:np),pp(np+1:np*2),up,vp, nx,ny,dx,dy,nu);
[sol,flg,res,iter] = bicgstab(afun,rhs, 1.0e-5, 1000);
disp(['solver: flag=',int2str(flg),'; res=',num2str(res), ';iter=',int2str(iter)]);

fp = sol(1:np);
gp = sol(np+1:np*2);

if (1)
    figure;
	plot(xp,yp,'.b');
	hold on;
    quiver(xp,yp,fp,gp);
	hold off;
    axis equal;
    axis([0 Lx 0 Ly]);
end

%
%
fs = Emat * fp;
gs = Emat * gp;
fs = reshape(fs, nx,ny);
gs = reshape(gs, nx,ny);

[us,vs,ps] = FastSolve2d(kx,ky,fs,gs, nx,ny,dx,dy,nu);

uper = FillPeriodic(us,nx,ny);
vper = FillPeriodic(vs,nx,ny);
pper = FillPeriodic(ps,nx,ny);


if 1
    figure;
    contourf(xper,yper,pper);
    
	hold on;
	% velocity
    % quiver(xs,ys,us,vs);
	%
	plot(xp,yp,'.k');
	
    hold off;
    axis equal;
    axis([0 Lx 0 Ly]);
end

if 1
	figure;
	upm = mean(up);
	vpm = mean(vp);
	quiver(xper,yper, uper-upm,vper-vpm);
	xstart = Lx*0.5 .* ones(ny,1);
	ystart = linspace(dy/2,Ly-dy/2,ny)';
	% streamline(stream2(xper',yper', uper'-upm,vper'-vpm, xstart,ystart));
	
	hold on;
	plot(xp,yp,'.k');
	hold off;
	
    axis equal;
    axis([0 Lx 0 Ly]);
end

% uana = cos(xs).*sin(ys);
% vana = -sin(xs).*cos(ys);
% pana = -0.25 * (cos(2*xs) + cos(2*ys));
