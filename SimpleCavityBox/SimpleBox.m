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

% ## SimpleBox

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-02

clc;
% clf;
clear all;

simple_globals;

% initialization
Lx = 1;
Ly = 1;
% ncell = 40;
ncell = 128;
nx = ncell;
ny = ncell;
dx = Lx / nx;
dy = Ly / ny;

relax = 0.8;
urelax = relax;
vrelax = relax;
prelax = 1 - relax;

ULid = 1;
rho = 1;
Re = 3200;
visc = rho*Lx*ULid / Re;

% dt = 0.025 * Re * dx^2;
% dt = 0.001;
dt = Inf;

% storage
umac = zeros(nx+1,ny+2);
vmac = zeros(nx+2,ny+1);
ustar = zeros(nx+1,ny+2);
vstar = zeros(nx+2,ny+1);
p = zeros(nx+2,ny+2);
pstar = zeros(nx+2,ny+2);
pdash = zeros(nx+2,ny+2);

umac = apply_umac_bc(umac);
vmac = apply_vmac_bc(vmac);

uold = umac;
vold = vmac;
pold = p;

% matrix coefficients
% Aw = zeros(nx,ny);
% As = zeros(nx,ny);
% An = zeros(nx,ny);
% Ae = zeros(nx,ny);
% Ap = zeros(nx,ny);
% bp = zeros(nx,ny);

%
AUp = zeros(nx+1,ny+2);
AVp = zeros(nx+2,ny+1);

% computation begins
disp(['Steady Lid-driven cavity, Re=', num2str(Re)]);

max_steps = 5000;
for istep = 1:max_steps
    
    uold = umac;
    vold = vmac;
    pold = pstar;
    
    ustar = solve_ustar(umac,vmac,uold,vold,pstar);
    vstar = solve_vstar(umac,vmac,uold,vold,pstar);
    
    [umac,vmac,pstar] = solve_pressure(ustar,vstar,pstar);
    
    % compute intermeidate velocity
    
    
    if (mod(istep,10) == 0)
        
        ucorr = umac - uold;
        vcorr = vmac - vold;
        % pcorr = pstar - pold;
        
        tol_abs = 1e-4 * ULid;
        corr = max([norm(ucorr(:),inf), norm(vcorr(:),inf)]);
        
        disp(['step=', int2str(istep), ...
            ', corr=', num2str(corr)]);
        
        if (corr < tol_abs)
            break;
        end
    end
end

xcs = linspace(dx/2,Lx-dx/2,nx);
ycs = linspace(dy/2,Ly-dy/2,ny);

ucell = 0.5 * (umac(1:nx,2:ny+1) + umac(2:nx+1,2:ny+1));
vcell = 0.5 * (vmac(2:nx+1,1:ny) + vmac(2:nx+1,2:ny+1));
pcell = pstar(2:nx+1,2:ny+1);

figure;
quiver(xcs, ycs, ucell', vcell');
title('Cell velocity');
% xlabel('x');
% ylabel('y');
axis equal;
axis([0 Lx 0 Ly]);

figure;
contourf(xcs,ycs,pcell',20);
axis equal;
axis([0 Lx 0 Ly]);
title("Pressure");

psi = easy_streamfunc(xcs,ycs,nx,ny,dx,dy,ucell,vcell);
figure;
contour(xcs, ycs, psi', 80);
title('stream-function');
axis equal;
axis([0 Lx 0 Ly]);





