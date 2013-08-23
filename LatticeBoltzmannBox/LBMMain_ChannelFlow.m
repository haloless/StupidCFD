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

% ## LBMMain_ChannelFlow

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-13

clc;
clear all;

refine = 1;
Lx = 400 * refine;
Ly = 100 * refine;
ObjCx = Lx/5;
% ObjCy = Ly/2 + 3*refine;
ObjCy = Ly/2;
ObjD = Ly/5;
ObjR = ObjD/2;

nx = Lx;
ny = Ly;
hx = Lx / nx;
hy = Ly / ny;
hh = min([hx hy]);
qc = 1;
dt = hh / qc;
cs2 = 1/3 * qc^2;

% determine parameters
Re = 100;
UMax = 0.1;
% % relaxation parameter, must < 2
% omega = 1.7;

if (exist('omega'))
    tau = 1 / omega;
    nu = (tau-0.5) * cs2 * dt;
    UMax = Re*nu/ObjD;
elseif (exist('UMax'))
    nu = UMax*ObjD / Re;
    tau = 1/cs2*nu/dt + 0.5;
    omega = 1/tau;
else
    error('Input parameters not given');
end

% compute some parameters
Re_h = UMax * hh / nu;
Cfl = dt*UMax / hh;
Dfl = dt*nu / hh^2;

% D2Q9
% [qwgt qex qey qord qopp] = LBMD2Q9Model();
qwgt = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36];
qex  = [0, 1, 0, -1, 0, 1, -1, -1, 1];
qey  = [0, 0, 1, 0, -1, 1, 1, -1, -1];
qord = [1, 2, 3, 4, 5,  6, 7, 8,  9];
qopp = [1, 4, 5, 2, 3,  8, 9, 6,  7];

[X,Y] = ndgrid(linspace(0+hx/2,Lx-hx/2,nx),linspace(0+hy/2,Ly-hy/2,ny));
% solid region
ebmask = zeros(nx,ny);
ebmask = (X-ObjCx).^2+(Y-ObjCy).^2 <= ObjR^2;
ebmask(:,1) = 1;
ebmask(:,ny) = 1;
ebRegion = find(ebmask);
% BC region
inletIdx = sub2ind([nx,ny], 1*ones(1,ny-2), 2:ny-1);
outletIdx = sub2ind([nx,ny], nx*ones(1,ny-2), 2:ny-1);

% initial condition with Poiseuille flow
L = Ly - 2*hy;
YPhys = Y - 1*hy;
ux = 4 * UMax/L^2 * (YPhys.*L - YPhys.^2);
% ux = UMax*ones(nx,ny);
UIn = ux(1,2:ny-1);
uy = zeros(nx,ny);
rho = ones(nx,ny);
% reshape to column vector
ux = reshape(ux, nx*ny,1);
uy = reshape(uy, nx*ny,1);
rho = reshape(rho, nx*ny,1);
% distrib. function
fs = zeros(nx*ny,9);
for i = 1:9
    eu = qex(i)*ux + qey(i)*uy;
    fs(:,i) = qwgt(i) * rho .* (1 + 3*eu + 9/2*eu.^2 - 3/2*(ux.^2+uy.^2));
end

max_step = 20000;
max_time = max_step * dt;
time = 0;
step = 0;
while (step<max_step && time<max_time)
    time = time + dt;
    step = step + 1;
    
    % macroscopic variables
    rho = sum(fs')';
    ux = (fs * qex') ./ rho;
    uy = (fs * qey') ./ rho;
    
    % macroscopic BC
    % inlet, Poiseuille
    % ux(inletIdx) = 4*UMax/L^2 * (YPhys(inletIdx)*L - YPhys(inletIdx).^2);
    ux(inletIdx) = UIn;
    uy(inletIdx) = 0;
    rho(inletIdx) = 1 ./ (1-ux(inletIdx)) .* ...
        (sum(fs(inletIdx,[1,3,5])')' + 2*sum(fs(inletIdx,[4,7,8])')');
    % outlet, constant pressure
    rho(outletIdx) = 1;
    ux(outletIdx) = -1 + 1 ./ rho(outletIdx) .* ...
        (sum(fs(outletIdx,[1,3,5])')' + 2*sum(fs(outletIdx,[2,6,9])')');
    uy(outletIdx) = 0;
    
    % microscopic BC
    % inlet, Zou-He BC
    fs(inletIdx,2) = fs(inletIdx,4) + 2/3*rho(inletIdx).*ux(inletIdx);
    fs(inletIdx,6) = fs(inletIdx,8) + 1/2*(fs(inletIdx,5)-fs(inletIdx,3)) ...
        + 1/2*rho(inletIdx).*uy(inletIdx) + 1/6*rho(inletIdx).*ux(inletIdx);
    fs(inletIdx,9) = fs(inletIdx,7) + 1/2*(fs(inletIdx,3)-fs(inletIdx,5)) ...
        - 1/2*rho(inletIdx).*uy(inletIdx) + 1/6*rho(inletIdx).*ux(inletIdx);
    % outlet, Zou-He BC
    fs(outletIdx,4) = fs(outletIdx,2) - 2/3*rho(outletIdx).*ux(outletIdx);
    fs(outletIdx,8) = fs(outletIdx,6) + 1/2*(fs(outletIdx,3)-fs(outletIdx,5)) ...
        - 1/2*rho(outletIdx).*uy(outletIdx) - 1/6*rho(outletIdx).*ux(outletIdx);
    fs(outletIdx,7) = fs(outletIdx,9) + 1/2*(fs(outletIdx,5)-fs(outletIdx,3)) ...
        + 1/2*rho(outletIdx).*uy(outletIdx) - 1/6*rho(outletIdx).*ux(outletIdx);
    
    % collision
    for i = 1:9
        eu = qex(i)*ux + qey(i)*uy;
        feq(:,i) = qwgt(i) * rho .* (1 + 3*eu + 9/2*eu.^2 - 3/2*(ux.^2+uy.^2));
        fout(:,i) = fs(:,i) - omega*(fs(:,i)-feq(:,i));
    end
    
    % bounce-back
    for i = 1:9
        fout(ebRegion,i) = fs(ebRegion,qopp(i));
    end
    
    % streaming
    for i = 1:9
        fs(:,i) = LBMStream2d(fout(:,i), nx,ny, qex(i),qey(i));
    end
    
    
    if (mod(step,50) == 0)
        prompt = ['step=',int2str(step)];
        disp(prompt);
        
        subnrow = 2;
        subncol = 1;
        
        velx = reshape(ux, nx,ny);
        vely = reshape(uy, nx,ny);
        
        subplot(subnrow,subncol,1);
        velmag = sqrt(velx.^2+vely.^2);
        velmag(ebRegion) = NaN;
        contourf(X',Y',velmag', 32);
        colorbar; 
        % shading flat;
        axis equal;
        title([prompt]);
        
        subplot(subnrow,subncol,2);
        pres = reshape(cs2*rho, nx,ny);
        pres(ebRegion) = NaN;
        contourf(X',Y',pres', 32);
        colorbar;
        axis equal;
        
        drawnow;
    end
end









