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

% restart = 1;
restart = 0;

if (~restart)
clc;
clear all;

refine = 1;
ObjConf = 0; % channel
% ObjConf = 1; % squares
% ObjConf = 2; % regular cylinder array
% W = 64;
% W = 50;
switch ObjConf
case {0, 1}
    W = 60;
    Lx = W * 10 * refine;
    Ly = W * 2 * refine;
    Re = 1;
    Bn = 10;
case {2}
    % W = 240;
    W = 180;
    % W = 120;
    Lx = W * 2 * refine;
    Ly = W * 1 * refine;
    Re = 1;
    Bn = 20;
otherwise
    error('Invalid Object Config.');
end

nx = Lx + 1;
ny = Ly + 1;
hx = Lx / (nx-1);
hy = Ly / (ny-1);
hh = min([hx hy]);
qc = 1;
dt = hh / qc;
cs = 1/sqrt(3) * qc;
cs2 = 1/3 * qc^2;


% determine parameters

%
rho0 = 1;
% UMax = 0.0002;
UMax = 0.00025;
etap = rho0 * UMax * (2*W) / Re;
tau0 = Bn * etap * UMax / W;
% constitutive parameters
% mreg = 1.0e10;
% mreg = 1.0e8;
% mreg = 5.0e6;
mreg = 4.0e6;
% mreg = 1.0e6;
eta0 = mreg*tau0 + etap;
% dotgamma0 = 1e-12; % min cutoff of shear rate
% eta0 = -tau0/dotgamma0*expm1(-mreg*dotgamma0) + etap;
dotgamma0 = 0;
lambda0 = (eta0/rho0) / (cs2*dt) + 0.5;


% if (exist('UMax'))
    % nu = UMax*ObjD / Re;
    % tau = 1/cs2*nu/dt + 0.5;
    % omega = 1/tau;
% else
    % error('Input parameters not given');
% end

% compute some parameters
% Re_h = UMax * hh / nu;
Cfl = dt*UMax / hh;
% Dfl = dt*nu / hh^2;

% D2Q9
% [qwgt qex qey qord qopp] = LBMD2Q9Model();
qwgt = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36];
qex  = [0, 1, 0, -1, 0, 1, -1, -1, 1];
qey  = [0, 0, 1, 0, -1, 1, 1, -1, -1];
qord = [1, 2, 3, 4, 5,  6, 7, 8,  9];
qopp = [1, 4, 5, 2, 3,  8, 9, 6,  7];

% [X,Y] = ndgrid(linspace(0+hx/2,Lx-hx/2,nx),linspace(0+hy/2,Ly-hy/2,ny));
xgrids = linspace(0,Lx,nx);
ygrids = linspace(0,Ly,ny);
[X,Y] = ndgrid(xgrids,ygrids);
% solid region
ebmask = zeros(nx,ny);
switch (ObjConf)
case {0} % nothing
    ebmask(:) = 0;
case {1} %
    ObjD = W * 2/3;
    ObjX = ObjD * 5/2; ObjY = Ly / 2;
    ebmask(abs(X-ObjX)<ObjD/2 & abs(Y-ObjY)<ObjD/2) = 1;
    ObjX = ObjD * 15/2; ObjY = Ly / 2;
    ebmask(abs(X-ObjX)<ObjD/2 & abs(Y-ObjY)<ObjD/2) = 1;
    ObjX = ObjD * 25/2; ObjY = Ly / 2;
    ebmask(abs(X-ObjX)<ObjD/2 & abs(Y-ObjY)<ObjD/2) = 1;
case {2} %
    ObjD = W / 8;
    ObjNx = 4; ObjNy = 3;
    for j = 1:ObjNy
    for i = 1:ObjNx
        ObjX = 1/2 * (Lx - ObjD*(ObjNx-1)*2) + (i-1)*ObjD*2;
        ObjY = 1/2 * (Ly - ObjD*(ObjNy-1)*2) + (j-1)*ObjD*2;
        ebmask((X-ObjX).^2+(Y-ObjY).^2 <= (ObjD^2)/4) = 1;
    end
    end
    % ebmask = (X-ObjCx).^2+(Y-ObjCy).^2 <= ObjR^2;
otherwise
    error('Invalid object configuration');
end
% channel wall
ebmask(:,1) = 1;
ebmask(:,ny) = 1;
ebRegion = find(ebmask);
% BC region
inletIdx = sub2ind([nx,ny], 1*ones(1,ny-2), 2:ny-1);
outletIdx = sub2ind([nx,ny], nx*ones(1,ny-2), 2:ny-1);

% initial condition with Poiseuille flow
if (0)
    L = Ly - 2*hy;
    YPhys = Y - 1*hy;
    ux = 4 * UMax/L^2 * (YPhys.*L - YPhys.^2);
else
    % ux = UMax*ones(nx,ny);
    ux = zeros(nx,ny);
    ux(inletIdx) = UMax;
end
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
feq = zeros(nx*ny,9);
fneq = zeros(nx*ny,9);

% non-Newton
lambda = zeros(nx*ny,1);
lambda(:) = lambda0;
%
dotgamma = zeros(nx*ny,1);
eta = zeros(nx*ny,1);
DMat = zeros(nx*ny,2*2);
%

time = 0;
step = 0;
end % not a restart

% max_step = 100000;
max_step = 50000;
max_time = max_step * dt;

while (step<max_step && time<max_time)
    time = time + dt;
    step = step + 1;
    
    % macroscopic variables
    rho = sum(fs,2);
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
    
    % equilibrium and non-equilibrium
    for i = 1:9
        eu = qex(i)*ux + qey(i)*uy;
        feq(:,i) = qwgt(i) * rho .* (1 + 3*eu + 9/2*eu.^2 - 3/2*(ux.^2+uy.^2));
        fneq(:,i) = fs(:,i) - feq(:,i);
        % fout(:,i) = fs(:,i) - omega*fneq(:,i);
    end
    % rate of deformation tensor
    DMat(:,1) = sum(fneq * diag(qex.*qex), 2);
    DMat(:,2) = sum(fneq * diag(qey.*qex), 2);
    DMat(:,3) = sum(fneq * diag(qex.*qey), 2);
    DMat(:,4) = sum(fneq * diag(qey.*qey), 2);
    if (1) % explicit
        % for i = 1:4
            % DMat(:,i) = -1/(2*cs2*dt) * DMat(:,i) ./ (rho.*lambda);
        % end
        %
        % dotgamma = sqrt(2*sum(DMat.^2,2));
        dotgamma = 2 * 1/(2*cs2*dt) ./ (rho.*lambda) .* sqrt(sum(DMat.^2,2));
        % non-Newtonian viscosity
        % Papanastasiou (or modified Bingham)
        if (1) % take care of zero shear
            % mask = (dotgamma ~= 0);
            mask = (dotgamma > dotgamma0);
            eta(mask) = -tau0 ./ dotgamma(mask) .* expm1(-mreg*dotgamma(mask)) + etap;
            eta(~mask) = eta0;
        else
            eta = -tau0 ./ dotgamma .* expm1(-mreg*dotgamma) + etap;
        end
        % relaxation time
        lambda = 1/(cs2*dt) * (eta ./ rho) + 0.5;
    else % implicit, Papanastasiou
        hatgamma = sqrt(sum(DMat.^2,2));
        mask0 = (hatgamma == 0);
        % mask0 = (hatgamma < 1e-12);
        hatgamma(hatgamma<1e-12) = 1e-12;
        res = zeros(size(lambda));
        jac = zeros(size(lambda));
        del = zeros(size(lambda));
        % initialize
        a = 1/(cs2*dt) ./ rho;
        % lambda(~mask0) = 0.5 + 1000*etap*a(~mask0);
        lambda(~mask0) = lambda0;
        % lambda(~mask0) = 0.5 + etap*a(~mask0);
        maxiter = 25;
        mask = mask0;
        for iter = 1:maxiter
            mask1 = (~mask);
            % mask1 = (~mask0);
            b = mreg * a ./ lambda;
            E1 = expm1(-hatgamma .* b);
            tmp1 = tau0 ./ hatgamma(mask1);
            %
            res(mask1) = lambda(mask1) ...
            + tmp1 .* lambda(mask1) .* E1(mask1) ...
            - 0.5 - etap*a(mask1);
            jac(mask1) = 1 + (tmp1 + tau0*b(mask1)) .* E1(mask1) + tau0*b(mask1);
            %
            del(mask1) = -res(mask1) ./ jac(mask1);
            lambda(mask1) = lambda(mask1) + del(mask1);
            
            res(:) = 0;
            res(mask1) = lambda(mask1) ...
            + tmp1 .* lambda(mask1) .* E1(mask1) ...
            - 0.5 - etap*a(mask1);
            mask = mask | (abs(res)<1e-10);
        end
        if (true && any(~mask))
            % find(~mask)
            length(find(~mask))
            warning('not fully converged')
        end
        lambda((~mask) & (lambda<0.5)) = 0.5;
        lambda(mask0) = lambda0;
        eta(mask0) = eta0;
        eta(~mask0) = (cs2*dt) * (lambda(~mask0)-0.5) .* rho(~mask0);
    end
    % collision
    omega = 1 ./ lambda;
    for i = 1:9
        fout(:,i) = fs(:,i) - omega .* fneq(:,i);
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
        prompt = ['step=',int2str(step), ';time=',num2str(time)];
        disp(prompt);
        
        subnrow = 3;
        subncol = 1;
        
        velx = reshape(ux, nx,ny);
        vely = reshape(uy, nx,ny);
        
        subplot(subnrow,subncol,1);
        velmag = sqrt(velx.^2+vely.^2);
        velmag(ebRegion) = NaN;
        % contourf(X',Y',velmag', 16);
        imagesc(xgrids,ygrids,velmag');
        colorbar; 
        % shading flat;
        axis equal;
        title([prompt]);
        
        subplot(subnrow,subncol,2);
        visc = reshape(eta./etap, nx,ny);
        if (ObjConf == 0)
            visc = log(visc);
        end
        visc(ebRegion) = NaN;
        % imagesc(xgrids,ygrids,visc', [1,10]);
        hpcolor = pcolor(X,Y,visc); caxis([1,10]);
        set(hpcolor,'EdgeColor','None');
        colorbar;
        axis equal;
        
        % pres = reshape(cs2*rho, nx,ny);
        % pres(ebRegion) = NaN;
        % contourf(X',Y',pres', 16);
        
        subplot(subnrow,subncol,3);
        DMat(:,1) = sum(fneq * diag(qex.*qex), 2);
        DMat(:,2) = sum(fneq * diag(qey.*qex), 2);
        DMat(:,3) = sum(fneq * diag(qex.*qey), 2);
        DMat(:,4) = sum(fneq * diag(qey.*qey), 2);
        for i = 1:4
            DMat(:,i) = -1/(2*cs2*dt) * DMat(:,i) ./ (rho.*lambda);
        end
        stress = 2 * eta .* sqrt(sum(DMat.^2,2));
        yield = zeros(size(eta));
        yield(:) = 2; % yielded
        yield(stress<=tau0) = 1; % un-yielded
        yield(ebRegion) = 0;
        yield = reshape(yield, nx,ny);
        imagesc(xgrids,ygrids,yield');
        colorbar;
        axis equal;
        
        drawnow;
    end
end









