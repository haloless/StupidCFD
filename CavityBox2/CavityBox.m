% ## Copyright (C) 2013 homu
% ## 
% ## This program is free software; you can redistribute it and/or modify
% ## it under the terms of the GNU General Public License as published by
% ## the Free Software Foundation; either version 3 of the License, or
% ## (at your option) any later version.
% ## 
% ## This program is distributed in the hope that it will be useful,
% ## but WITHOUT ANY WARRANTY; without even the implied warranty of
% ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% ## GNU General Public License for more details.
% ## 
% ## You should have received a copy of the GNU General Public License
% ## along with Octave; see the file COPYING.  If not, see
% ## <http://www.gnu.org/licenses/>.

% ## main

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-06-26

% clc
% clf
clear all

% ncell = 10;
% ncell = 15;
% ncell = 20;
% ncell = 30;
% ncell = 40;
% ncell = 50;
% ncell = 60;
% ncell = 80;
ncell = 100;

nx = ncell;
ny = ncell;
nxy = [nx, ny];
L = pi * 2;
% L = pi;
hx = L / nx;
hy = L / ny;
hh = min(hx, hy);
hxy = [hx, hy];

bc = struct; % dummy
bc.is_confined = 1;


% conditions
global ULid
ULid = 1; 
Re = 10.0;
nu = 1 / Re;
CFL = 0.25;
% dt = min([0.25*Re*hh^2, CFL*hh/ULid]);
dt = 0.50;
% dt = 0.40;
% dt = 0.20;
% dt = 0.10;
% dt = 0.08;
% dt = 0.05; % explicit stable
% dt = 0.04; 
% dt = 0.02;
% dt = 0.01;
% dt = 0.008;
% dt = 0.005;
% dt = 0.002;
% dt = 0.001;
% dt = 0.0005;
% dt = 0.0002;



% storages, forward staggered
p = zeros(nx+2, ny+2);
umac = zeros(nx+2,ny+2);
vmac = zeros(nx+2,ny+2);

for i = 1:nx
for j = 1:ny
	xx = (i-1)*hx;
	yy = (j-0.5)*hy;
	umac(i+1,j+1) = cos(xx) * sin(yy);
end
end
for i = 1:nx
for j = 1:ny
	xx = (i-0.5)*hx;
	yy = (j-1)*hy;
	vmac(i+1,j+1) = -sin(xx) * cos(yy);
end
end
[umac, vmac] = apply_mac_bc(umac,vmac, nx,ny, bc);


p_old = p;
umac_old = umac;
vmac_old = vmac;
uref = umac;
vref = vmac;

% intermediate values
ustar = umac;
vstar = vmac;
Hx = umac;
Hy = vmac;
Hx_old = Hx;
Hy_old = Hy;

% node position
cellXs = (-1:nx)*hx + hx/2;
cellYs = (-1:ny)*hy + hy/2;
edgeXs = (0:nx) * hx;
edgeYs = (0:ny) * hy;
% disp(cellXs);
% disp(cellYs);
% disp(edgeXs);
% disp(edgeYs);



Lap = mac_laplacian(nx,ny, hx,hy, bc);
perm = symamd(Lap);
RLap = chol(Lap(perm,perm));
% RLap = chol(Lap);
RLapt = RLap';

disp('MAC Laplacian OP built.');

visc_coef = 0.5;
Visc = mac_viscop(nx,ny,hx,hy, nu,dt*(1-visc_coef), bc);
disp('Visc OP built.');


% max_time = 0.01;
% max_time = 0.1;
max_time = 1.0;
% max_time = 5.0;
curr_time = 0;
% max_steps = 10000;
max_steps = round(max_time/dt);

rnorm = 0;

for istep = 1:max_steps
    curr_time = curr_time + dt;
    [umac, vmac] = apply_mac_bc(umac,vmac, nx,ny, bc);
    
	if 0
		[Hx, Hy] = mac_predictor(umac,vmac, nx,ny, hx,hy, dt, Re, bc);
		% if (istep == 1)
			% ustar = umac + dt * Hx;
			% vstar = vmac + dt * Hy;
		% else
			% ustar = umac + dt/2 * (3*Hx - Hx_old);
			% vstar = vmac + dt/2 * (3*Hy - Hy_old);
		% end
		ustar = umac + dt * Hx;
		vstar = vmac + dt * Hy;
    else
		[Hx, Hy, Lx, Ly] = mac_advect(umac,vmac, nx,ny, hx,hy, dt, Re, bc);
		utmp = umac + dt * Hx + visc_coef*dt*Lx;
		vtmp = vmac + dt * Hy + visc_coef*dt*Ly;
		
		utmp = reshape(utmp(2:nx+1,2:ny+1),[],1);
		vtmp = reshape(vtmp(2:nx+1,2:ny+1),[],1);
		ustar(2:nx+1,2:ny+1) = reshape(Visc \ utmp, nx,ny);
		vstar(2:nx+1,2:ny+1) = reshape(Visc \ vtmp, nx,ny);
	end
	
	[ustar,vstar] = apply_mac_bc(ustar,vstar, nx,ny, bc);
    
    % solve PPE
    rhs =  mac_rhs(ustar,vstar, nx,ny, hx,hy,dt);
	if 1
		% sol = Lap \ rhs;
		sol = rhs;
		sol(perm) = RLap \ (RLapt \ rhs(perm));
		% sol = RLap \ (RLapt \ rhs);
	else
		sol = bicgstab(Lap, rhs, 1.0e-8, 1000); 
	end
    
    phi = reshape(sol, nx,ny);
    % figure
    % contourf(cellXs(2:nx+1),cellYs(2:ny+1), phi, 20);
    % axis equal
    % axis([0 L 0 L]);
    % title('Phi')
    % phi(:) = 0;
    
    p(2:nx+1,2:ny+1) = phi;
    p = apply_pres_bc(p, nx, ny, bc);
    
    % correction
	irange = 2:nx+1;
	jrange = 2:ny+1;
	umac(irange,jrange) = ustar(irange,jrange) - dt/hx * (p(irange,jrange)-p(irange-1,jrange));
	vmac(irange,jrange) = vstar(irange,jrange) - dt/hy * (p(irange,jrange)-p(irange,jrange-1));
    % umac(2:nx,2:ny+1) = ustar(2:nx,2:ny+1) - ...
        % dt/hx * (p(3:nx+1,2:ny+1) - p(2:nx,2:ny+1));
    % vmac(2:nx+1,2:ny) = vstar(2:nx+1,2:ny) - ...
        % dt/hy * (p(2:nx+1,3:ny+1) - p(2:nx+1,2:ny));
    
    % for i = 2:nx
        % for j = 2:ny+1
            % umac(i,j) = ustar(i,j) - dt/hx * (p(i+1,j)-p(i,j));
        % end
    % end
    % for i = 2:nx+1
        % for j = 2:ny
            % vmac(i,j) = vstar(i,j) - dt/hy * (p(i,j+1)-p(i,j));
        % end
    % end
    [umac,vmac] = apply_mac_bc(umac, vmac, nx, ny, bc);

    
    if mod(istep,10) == 0 || istep==max_steps
        disp(['step = ', int2str(istep), ...
            '; time = ', num2str(curr_time)]);
        
        if mod(istep,100) == 0 || istep==max_steps
            contourf(cellXs(2:nx+1),cellYs(2:ny+1), phi', 20);
            axis equal;
            axis([0 L 0 L]);
			colorbar;
            drawnow;
        end
        
    end
    
    % Hx_old = ustar - umac_old;
    % Hy_old = vstar - vmac_old;
    Hx_old = Hx;
    Hy_old = Hy;
    umac_old = umac;
    vmac_old = vmac;
	
	% return
	
	if 0
		% calculate max error of L2 error
		cana = exp(-2*nu*curr_time);
		uana = uref .* cana;
		vana = vref .* cana;
		uerr = umac(2:nx+1,2:ny+1) - uana(2:nx+1,2:ny+1);
		% verr = vmac(2:nx+1,2:ny+1) - vana(2:nx+1,2:ny+1);
		err = norm(uerr(:));
		% err = max(abs(uerr(:)));
		rnorm = max(rnorm, err);
	end
end

cana = exp(-2*nu*curr_time);
uana = uref .* cana;
vana = vref .* cana;
uerr = umac(2:nx+1,2:ny+1) - uana(2:nx+1,2:ny+1);
verr = vmac(2:nx+1,2:ny+1) - vana(2:nx+1,2:ny+1);
%
load('usol.mat');
uerr = umac(2:nx+1,2:ny+1) - usol(2:nx+1,2:ny+1);

% rnorm = sqrt(sum(uerr(:).^2 + verr(:).^2) / (nx*ny));
% rnorm = norm(uerr(:)) + norm(verr(:));
rnorm = max(abs(uerr(:)));
dt
rnorm
% rnorm / ncell


return













% post processing
% compute cell velocity, pressure, position
ucell = 0.5 * (umac(1:nx,2:ny+1) + umac(2:nx+1,2:ny+1));
vcell = 0.5 * (vmac(2:nx+1,1:ny) + vmac(2:nx+1,2:ny+1));
pcell = p(2:nx+1,2:ny+1);
xs = cellXs(2:nx+1);
ys = cellYs(2:ny+1);

% plot velocity vectors
figure;
quiver(xs, ys, ucell', vcell');
title('Cell velocity');
xlabel('x');
ylabel('y');
axis equal;
axis([0 L 0 L]);

% streamline
psi = easy_streamfunc(xs,ys,nx,ny,hx,hy,ucell,vcell);
figure;
contour(xs, ys, psi', 100);
title('stream-function');
axis equal;
axis([0 L 0 L]);

% yet another streamline
% Npsi = (nx-1)*(ny-1);
% % Lpsi = spalloc(Npsi,Npsi,Npsi*5);
% Lpsi = kron(speye(ny-1),K1(nx-1,hx,2))+kron(K1(ny-1,hy,2),speye(nx-1));
% rhs = zeros(nx-1,ny-1);
% rhs = rhs - 1/hy * (umac(2:nx,3:ny+1) - umac(2:nx,2:ny));
% rhs = rhs + 1/hx * (vmac(3:nx+1,2:ny) - vmac(2:nx,2:ny));
% rhs = reshape(rhs,Npsi);
% psi2 = Lpsi \ rhs;
% psi2 = reshape(psi2,nx-1,ny-1);
psi2 = poormans_streamfunc(nx,ny,hx,hy,umac,vmac);
figure;
contour(0.5*(xs(1:end-1)+xs(2:end)),0.5*(ys(1:end-1)+ys(2:end)),psi2',100);
title('stream-function2');
axis equal;
axis([0 L 0 L]);

% plot pressure contour
figure;
[C,h] = contourf(xs, ys, p(2:nx+1,2:ny+1)', 20);
axis equal
axis([0 L 0 L]);
title('Pressure')
% clabel(C,h)

