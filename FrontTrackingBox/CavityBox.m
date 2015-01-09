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

clc
% clf
clear all

ncell = 100;
nx = ncell;
ny = ncell;
nxy = [nx, ny];
L = 1;
hx = L / nx;
hy = L / ny;
hh = min(hx, hy);
hxy = [hx, hy];

% storages, forward staggered
p = zeros(nx+2, ny+2);
umac = zeros(nx+1,ny+2);
vmac = zeros(nx+2,ny+1);
p_old = p;
umac_old = umac;
vmac_old = vmac;

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

% conditions
global ULid
ULid = 1; % lid velocity
Re = 3200;
CFL = 0.25;
dt = min([0.25*Re*hh^2, CFL*hh/ULid]);

bc = struct; % dummy
bc.is_confined = 1;

Lap = mac_laplacian(nx,ny, hx,hy, bc);
perm = symamd(Lap);
RLap = chol(Lap(perm,perm));
% RLap = chol(Lap);
RLapt = RLap';

disp('MAC Laplacian OP built.');

curr_time = 0;
max_steps = 20000;
for istep = 1:max_steps
    curr_time = curr_time + dt;
    % [umac, vmac] = apply_mac_bc(umac,vmac, nx,ny, bc);
    
    [Hx, Hy] = mac_predictor(umac,vmac, nx,ny, hx,hy, dt, Re, bc);
    if (istep == 1)
    % if(istep > 0)
        ustar = umac + dt * Hx;
        vstar = vmac + dt * Hy;
    else
        ustar = umac + dt/2 * (3*Hx - Hx_old);
        vstar = vmac + dt/2 * (3*Hy - Hy_old);
    end
    [ustar,vstar] = apply_mac_bc(ustar,vstar, nx,ny, bc);
    
    % solve PPE
    rhs =  mac_rhs(ustar,vstar, nx,ny, hx,hy,dt);
    % sol = Lap \ rhs;
    sol = rhs;
    sol(perm) = RLap \ (RLapt \ rhs(perm));
    % sol = RLap \ (RLapt \ rhs);
    
    phi = reshape(sol, nx,ny);
    % figure
    % contourf(cellXs(2:nx+1),cellYs(2:ny+1), phi, 20);
    % axis equal
    % axis([0 L 0 L]);
    % title('Phi')
    
    
    p(2:nx+1,2:ny+1) = phi;
    p = apply_pres_bc(p, nx, ny, bc);
    
    % correction
    umac(2:nx,2:ny+1) = ustar(2:nx,2:ny+1) - ...
        dt/hx * (p(3:nx+1,2:ny+1) - p(2:nx,2:ny+1));
    vmac(2:nx+1,2:ny) = vstar(2:nx+1,2:ny) - ...
        dt/hy * (p(2:nx+1,3:ny+1) - p(2:nx+1,2:ny));
    
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

    
    if mod(istep,10) == 0
        disp(['step = ', int2str(istep), ...
            '; time = ', num2str(curr_time)]);
        
        if mod(istep,100) == 0
            contourf(cellXs(2:nx+1),cellYs(2:ny+1), phi', 20);
            axis equal;
            axis([0 L 0 L]);
            drawnow;
        end
        
        eps_abs = 1.0e-5;
        umac_diff = umac - umac_old;
        err = norm(umac_diff(:),Inf);
        disp(['err = ', num2str(err)]);
        if(err <= eps_abs)
            break;
        end
    end
    
    % Hx_old = ustar - umac_old;
    % Hy_old = vstar - vmac_old;
    Hx_old = Hx;
    Hy_old = Hy;
    umac_old = umac;
    vmac_old = vmac;
end

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

