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

% ## LBM_RayleighBenardInstability

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-23

% see http://onlinelibrary.wiley.com/doi/10.1002/fld.337/abstract

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

Ly = 1;
Lx = Ly*2;
refine = 1;
% nx = refine*128 + 2;
nx = refine*128;
ny = refine*64 + 2;
dh = Ly / (ny-2);
%
% xs = linspace(-dh/2,Lx+dh/2,nx);
xs = linspace(dh/2,Lx-dh/2,nx);
ys = linspace(-dh/2,Ly+dh/2,ny);
[X,Y] = ndgrid(xs,ys);
%
idx_ylo = sub2ind([nx,ny], 1:nx, repmat(1,1,nx));
idx_yhi = sub2ind([nx,ny], 1:nx, repmat(ny,1,nx));

%
Pr = 1;
Ra = 20000; % Rayleigh No.
gr = 0.001;
buoyancy = [0, -gr];
%
THot = 1; % heating on bottom wall
TCold = 0; % cooling on top wall
T0 = 1/2 * (THot + TCold);

% parameters
dt = sqrt(gr*dh);
% kinematic viscosity in lattice unit
nu = sqrt(Pr/Ra) * dt/dh^2;
% thermal diffusion in lattice unit
kappa = 1/sqrt(Pr*Ra) * dt/dh^2;

% relaxation
tauNS = 3*nu + 0.5;
tauT = 3*kappa + 0.5;
omegaNS = 1/tauNS;
omegaT = 1/tauT;

% storage
rho = ones(nx*ny,1);
u = zeros(nx*ny,1);
v = zeros(nx*ny,1);
%
T = TCold * ones(nx,ny);
T(:,1) = THot;
% perturb the temperature field
T(nx/2,2) = THot * 1.1;
T = reshape(T, nx*ny,1);

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
    gBous = -1/(THot-TCold) * (rho .* (T-T0)) * buoyancy;
    for i = 1:9
        eu = exNS(i)*u + eyNS(i)*v;
        fEq(:,i) = wNS(i) * rho .* (1 + 3*eu + 9/2*eu.^2 - 3/2*vel2);
        fExt(:,i) = 3*wNS(i) * (gBous * [exNS(i); eyNS(i)]);
        % fStar(:,i) = fs(:,i) + omegaNS*(fEq(:,i)-fs(:,i)) + fExt(:,i);
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
    % y-low
    fs(idx_ylo,3) = fs(idx_ylo,5);
    fs(idx_ylo,6) = fs(idx_ylo,8);
    fs(idx_ylo,7) = fs(idx_ylo,9);
    % y-high
    fs(idx_yhi,5) = fs(idx_yhi,3);
    fs(idx_yhi,8) = fs(idx_yhi,6);
    fs(idx_yhi,9) = fs(idx_yhi,7);
    
    % BC
    % y-low
    ts(idx_ylo,3) = THot - sum(ts(idx_ylo,[1,2,4,5]),2);
    % y-high
    ts(idx_yhi,5) = TCold - sum(ts(idx_yhi,[1,2,3,4]),2);
    
    % macroscopic variables
    rho = sum(fs, 2);
    u = (fs * exNS') ./ rho;
    v = (fs * eyNS') ./ rho;
    T = sum(ts, 2);
    % BC
    u(idx_ylo) = 0;
    v(idx_ylo) = 0;
    u(idx_yhi) = 0;
    v(idx_yhi) = 0;
    
    if (mod(step,100)==0)
        velx = reshape(u, nx,ny);
        vely = reshape(v, nx,ny);
        temp = reshape(T, nx,ny);
        
        Nusselt = 1 + sum(v.*T) / (nx*kappa*(THot-TCold));
        
        prompt = ['step=',int2str(step), ...
            ';Nusselt=',num2str(Nusselt), ...
        ];
        disp(prompt);
        
        subplot(2,1,1);
        velmag = sqrt(velx.^2 + vely.^2);
        contourf(velmag', 32);
        colorbar; shading flat;
        title(prompt);
        
        hold on;
        quiver(velx',vely');
        hold off;
        
        
        subplot(2,1,2);
        contourf(temp', 32);
        colorbar; shading flat;
        title(prompt);
        
        drawnow;
    end
end






