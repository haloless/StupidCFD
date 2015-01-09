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


% global data
FTGlobals;

% 2D Cart.
% CoordSys = CoordSys_Cart;

% RZ coord.
CoordSys = CoordSys_RZ; 
ProbLo = [0.0; 0.0];
ProbHi = [1.0; 8.0];
ProbLen = ProbHi - ProbLo;


% nx = 4;
% ny = 32;

nx = 16;
ny = 128;
dx = ProbLen(1) / nx;
dy = ProbLen(2) / ny;

%
CellXs = linspace(ProbLo(1)-0.5*dx, ProbHi(1)+0.5*dx, nx+2)';
CellYs = linspace(ProbLo(2)-0.5*dy, ProbHi(2)+0.5*dy, ny+2)';
EdgeXs = linspace(ProbLo(1), ProbHi(1), nx+1)';
EdgeYs = linspace(ProbLo(2), ProbHi(2), ny+1)';
if CoordSys == CoordSys_Cart % Cart
    EdgeRs = ones(nx+1,1);
    CellRs = ones(nx+2,1);
elseif CoordSys == CoordSys_RZ % RZ
    EdgeRs = EdgeXs;
    CellRs = CellXs;
end


% storages, forward staggered
p = zeros(nx+2, ny+2);
umac = zeros(nx+1,ny+2);
vmac = zeros(nx+2,ny+1);
pold = p;
uold = umac;
vold = vmac;

% intermediate values
ustar = umac;
vstar = vmac;
Hx = umac;
Hy = vmac;
Hx_old = Hx;
Hy_old = Hy;


% conditions
global UIn
UIn = 1.0;
Re = 10;
CFL = 0.25;
dt = min([0.05*Re*dx^2, CFL*dx/UIn]);

bc = struct; % dummy
bc.is_confined = 1;

Lap = mac_laplacian(nx,ny, dx,dy, bc);
if (1)
perm = symamd(Lap);
RLap = chol(Lap(perm,perm));
% RLap = chol(Lap);
RLapt = RLap';
end

disp('MAC Laplacian OP built.');


vmac(2:nx+1,1:ny+1) = UIn;

curr_time = 0;
max_steps = 5000;
for istep = 1:max_steps
    curr_time = curr_time + dt;
    [umac, vmac] = apply_mac_bc(umac,vmac, nx,ny, bc);
    
    [Hx, Hy] = mac_predictor(umac,vmac, nx,ny, dx,dy, dt, Re, bc);
    ustar = umac + dt * Hx;
    vstar = vmac + dt * Hy;
    % if (istep == 1)
    % % if(istep > 0)
        % ustar = umac + dt * Hx;
        % vstar = vmac + dt * Hy;
    % else
        % ustar = umac + dt/2 * (3*Hx - Hx_old);
        % vstar = vmac + dt/2 * (3*Hy - Hy_old);
    % end
    [ustar,vstar] = apply_mac_bc(ustar,vstar, nx,ny, bc);
    
    % solve PPE
    rhs =  mac_rhs(ustar,vstar, nx,ny, dx,dy,dt);
    if (0)
    sol = Lap \ rhs;
    else
    sol = rhs;
    sol(perm) = RLap \ (RLapt \ rhs(perm));
    % sol = RLap \ (RLapt \ rhs);
    end
    
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
        dt/dx * (p(3:nx+1,2:ny+1) - p(2:nx,2:ny+1));
    vmac(2:nx+1,2:ny) = vstar(2:nx+1,2:ny) - ...
        dt/dy * (p(2:nx+1,3:ny+1) - p(2:nx+1,2:ny));
    
    
    if mod(istep,10) == 0
        prompt = ['step = ', int2str(istep), ...
            '; time = ', num2str(curr_time)];
        disp(prompt);
        
        if mod(istep,100) == 0
            xs = CellXs(2:nx+1);
            ys = CellYs(2:ny+1);
            
            % pressure
            subplot(1,3,1)
            contourf(xs,ys, phi', 20);
            axis equal;
            axis([ProbLo(1) ProbHi(1) ProbLo(2) ProbHi(2)]);
            title(prompt);
            
            % velocity vectors
            subplot(1,3,2);
            quiver(EdgeXs,EdgeYs, ...
                0.5 * (umac(1:nx+1,2:ny+2)+umac(1:nx+1,1:ny+1))', ...
                0.5 * (vmac(2:nx+2,1:ny+1)+vmac(1:nx+1,1:ny+1))');
            % quiver(xs, ys, ucell', vcell');
            % contourf(xs,ys, sqrt(ucell.^2+vcell.^2)', 10);
            title('Node velocity');
            xlabel('x');
            ylabel('y');
            axis equal;
            axis([ProbLo(1) ProbHi(1) ProbLo(2) ProbHi(2)]);
            
            
            % subplot(1,3,3);
            % contourf(xs,ys, sqrt(ucell.^2+vcell.^2)', 10);
            % title('Cell velocity');
            % xlabel('x');
            % ylabel('y');
            % axis equal;
            % axis([ProbLo(1) ProbHi(1) ProbLo(2) ProbHi(2)]);
            
            drawnow;
        end
        
        eps_abs = 1.0e-5;
        umac_diff = umac - uold;
        vmac_diff = vmac - vold;
        err = norm(vmac_diff(:),Inf);
        disp(['err = ', num2str(err)]);
        if(err <= eps_abs)
            break;
        end
    end
    
    % Hx_old = ustar - umac_old;
    % Hy_old = vstar - vmac_old;
    % Hx_old = Hx;
    % Hy_old = Hy;
    uold = umac;
    vold = vmac;
end

% post processing
% compute cell velocity, pressure, position
ucell = 0.5 * (umac(1:nx,2:ny+1) + umac(2:nx+1,2:ny+1));
vcell = 0.5 * (vmac(2:nx+1,1:ny) + vmac(2:nx+1,2:ny+1));
pcell = p(2:nx+1,2:ny+1);
xs = CellXs(2:nx+1);
ys = CellYs(2:ny+1);

if (0)
% plot velocity vectors
figure;
% quiver(xs, ys, ucell', vcell');
contourf(xs,ys, sqrt(ucell.^2+vcell.^2)', 10);
title('Cell velocity');
xlabel('x');
ylabel('y');
axis equal;
axis([0 ProbLen(1) 0 ProbLen(2)]);
end

if (0)
% plot pressure contour
figure;
[C,h] = contourf(xs, ys, p(2:nx+1,2:ny+1)', 20);
axis equal
axis([0 ProbLen(1) 0 ProbLen(2)]);
title('Pressure')
% clabel(C,h)
end

if CoordSys == CoordSys_Cart
% analytical solution (2D)
xs = CellXs(2:nx+1);
vana = 6*UIn/4*(1- xs.^2);
figure; plot(xs,vmac(2:nx+1,ny+1),'x', xs,vana);
legend('sim','ana')

elseif CoordSys == CoordSys_RZ
% analytical solution (axisymmetric RZ)
rs = CellRs(2:nx+1);
vana = 2*UIn * (1-rs.^2);
figure; plot(rs,vmac(2:nx+1,ny+1),'x', ...
rs,vmac(2:nx+1,96),'o', rs,vmac(2:nx+1,64),'s', ...
rs,vana);
legend('sim,1','sim,3/4','sim,1/2','ana');
end