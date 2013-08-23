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

% ## LBM_NaturalConvectionMain

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-23

clc; clear all;

% D2Q9 for NS equation
qNS = 9;
wNS = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36];
exNS  = [0, 1, 0, -1, 0, 1, -1, -1, 1];
eyNS  = [0, 0, 1, 0, -1, 1, 1, -1, -1];
% D2Q5 for thermal equation
qT = 5;
wT = [1/3, 1/6, 1/6, 1/6, 1/6];
exT = [0, 1, 0, -1, 0];
eyT = [0, 0, 1, 0, -1];

% L = 1;
% Lx = L;
% Ly = L;
% refine = 1;
% nx = refine*128 + 2;
% ny = refine*128 + 2;
refine = 1;
L = 128*refine;
Lx = L;
Ly = L;
nx = Lx + 2;
ny = Ly + 2;
dx = Lx / (nx-2);
dy = Ly / (ny-2);
dh = min([dx,dy]);
qc = 1;
dt = dh / qc;
cs2 = 1/3 * qc^2;
cs = sqrt(cs2);

%
xs = linspace(-dx/2,Lx+dx/2,nx);
ys = linspace(-dy/2,Ly+dy/2,ny);
[X,Y] = ndgrid(xs,ys);
%
idx_xlo = sub2ind([nx,ny], repmat(1,1,ny), 1:ny);
idx_xhi = sub2ind([nx,ny], repmat(nx,1,ny),1:ny);
idx_ylo = sub2ind([nx,ny], 1:nx, repmat(1,1,nx));
idx_yhi = sub2ind([nx,ny], 1:nx, repmat(ny,1,nx));

%
TH = 1; % heating on left wall
TC = 0; % cooling on right wall
T0 = 1/2 * (TH + TC);
DT = TH - TC;
%
gr = 0.1;
buoyancy = [0, -gr];
%
Pr = 0.71; % Pr of air
Ra = 1e5; % Rayleigh 
%
omegaNS = 1.2;
tauNS = 1/omegaNS;
nu = (tauNS-0.5) * cs2 * dt;
%
kappa = nu / Pr;
tauT = kappa/(cs2*dt) + 0.5;
omegaT = 1/tauT;
%
beta = Ra*nu*kappa / (gr*DT*L^3);


% storage
rho = ones(nx*ny,1);
u = zeros(nx*ny,1);
v = zeros(nx*ny,1);
%
T = TC * ones(nx*ny,1);
T(idx_xlo) = TH;

%
fs = rho * wNS;
ts = T * wT;
%
fEq = zeros(nx*ny,9);
fExt = zeros(nx*ny,9);
fStar = zeros(nx*ny,9);
tEq = zeros(nx*ny,5);
tStar = zeros(nx*ny,5);

max_step = 40000;
max_time = max_step*dt;
step = 0;
time = 0;
while (step<max_step && time<max_time)
    step = step + 1;
    time = time + dt;
    
    % collision
    vel2 = u.^2 + v.^2;
    % Boussinesq
    gBous = -beta * (rho .* (T-T0)) * buoyancy;
    for i = 1:9
        eu = exNS(i)*u + eyNS(i)*v;
        fEq(:,i) = wNS(i) * rho .* (1 + 3*eu + 9/2*eu.^2 - 3/2*vel2);
        fExt(:,i) = 3*wNS(i) * (gBous * [exNS(i); eyNS(i)]);
    end
    fStar = (1-omegaNS)*fs + omegaNS*fEq + fExt;
    
    % collision
    for i = 1:5
        eu = exT(i)*u + eyT(i)*v;
        tEq(:,i) = wT(i) * T .* (1 + 3*eu);
    end
    tStar = (1-omegaT)*ts + omegaT*tEq;
    
    % streaming
    for i = 1:9
        fs(:,i) = LBMStream2d(fStar(:,i), nx,ny, exNS(i),eyNS(i));
    end
    for i = 1:5
        ts(:,i) = LBMStream2d(tStar(:,i), nx,ny, exT(i),eyT(i));
    end
    
    % BC
    % x-low: u=0, v=0, T=TH
    fs(idx_xlo,2) = fs(idx_xlo,4);
    fs(idx_xlo,6) = fs(idx_xlo,8);
    fs(idx_xlo,9) = fs(idx_xlo,7);
    ts(idx_xlo,2) = TH - sum(ts(idx_xlo,[1,3,4,5]),2);
    % x-high: u=0, v=0, T=TC
    fs(idx_xhi,4) = fs(idx_xhi,2);
    fs(idx_xhi,7) = fs(idx_xhi,9);
    fs(idx_xhi,8) = fs(idx_xhi,6);
    ts(idx_xhi,4) = TC - sum(ts(idx_xhi,[1,2,3,5]),2);
    % y-low: u=0, v=0, dT/dy=0
    fs(idx_ylo,3) = fs(idx_ylo,5);
    fs(idx_ylo,6) = fs(idx_ylo,8);
    fs(idx_ylo,7) = fs(idx_ylo,9);
    ts(idx_ylo,3) = ts(idx_ylo+nx,3);
    % y-high: u=0, v=0, dT/dy=0
    fs(idx_yhi,5) = fs(idx_yhi,3);
    fs(idx_yhi,8) = fs(idx_yhi,6);
    fs(idx_yhi,9) = fs(idx_yhi,7);
    ts(idx_yhi,5) = ts(idx_yhi-nx,5);
    
    % macroscopic variables
    rho = sum(fs, 2);
    u = (fs * exNS') ./ rho;
    v = (fs * eyNS') ./ rho;
    T = sum(ts, 2);
    % BC
    u(idx_xlo) = 0; v(idx_xlo) = 0; 
    u(idx_xhi) = 0; v(idx_xhi) = 0; 
    u(idx_ylo) = 0; v(idx_ylo) = 0;
    u(idx_yhi) = 0; v(idx_yhi) = 0;
    % T(idx_xlo) = TH;
    % T(idx_xhi) = TC;
    
    
    if (mod(step,100)==0)
        velx = reshape(u, nx,ny);
        vely = reshape(v, nx,ny);
        temp = reshape(T, nx,ny);
        
        % Nusselt = 1 + sum(v.*T) / (nx*kappa*(THot-TCold));
        
        prompt = ['step=',int2str(step), ...
        ];
        disp(prompt);
        
        subplot(1,2,1);
        % velmag = sqrt(velx.^2 + vely.^2);
        % contourf(velmag', 32);
        % colorbar; shading flat;
        psi = easy_streamfunc(xs,ys,nx,ny,dx,dy,velx,vely);
        contourf(psi', 32); colorbar;
        axis([0 Lx 0 Ly]); axis equal;
        title(prompt);
        
        subplot(1,2,2);
        contourf(temp', 32); colorbar;
        % shading flat;
        axis([0 Lx 0 Ly]); axis equal;
        title(prompt);
        
        drawnow;
    end
end





