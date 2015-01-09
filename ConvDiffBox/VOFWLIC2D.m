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

% ## VOFWLIC2D

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-10-03

clc; clear all;

Lx = 1.0;
Ly = 1.0;
ncell = 100;
nx = ncell;
ny = ncell;
dx = Lx / nx;
dy = Ly / ny;

cellxs = linspace(-dx/2,Lx+dx/2,nx+2);
cellys = linspace(-dy/2,Ly+dy/2,ny+2);
edgexs = linspace(-dx,Lx+dx,nx+3);
edgeys = linspace(-dy,Ly+dy,ny+3);
%
[X,Y] = ndgrid(cellxs,cellys);


f = zeros(nx+2,ny+2);
f0 = zeros(nx+2,ny+2);
umac = zeros(nx+3,ny+2);
vmac = zeros(nx+2,ny+3);

% initialize Zalesak Disk
% steps per single cycle 
cycleStep = 800;
omega = 2*pi / cycleStep;
% omega = 1.0;
% MAC velocity on cell edges
[xmac,ymac] = ndgrid(edgexs,cellys);
umac(:,:) = omega * (ymac - 0.5);
[xmac,ymac] = ndgrid(cellxs,edgeys);
vmac(:,:) = omega * (0.5 - xmac);
clear xmac ymac;
% VOF fraction
% TODO volume stack
f(:,:) = 1.0 * ((X-0.5).^2+(Y-0.75).^2<=0.17^2 & (abs(X-0.5)>0.03 | Y>0.85));
f0 = f;
f0sum = sum(f0(:));

if(1)
    figure;
    contourf(X',Y',f',10); axis equal;
    hold on;
    quiver(X',Y',0.5*(umac(1:nx+2,:)+umac(2:nx+3,:))',0.5*(vmac(:,1:ny+2)+vmac(:,2:ny+3))');
    hold off;
    xlabel('X'); ylabel('Y');
end

dt = 1.0;
% dt = dx;
time = 0.0;
maxStep = cycleStep;
% maxStep = 550;
for step = 1:maxStep
    time = time + dt;
    
    eps_small = 1e-10;
    f(f<eps_small) = 0;
    f(f>1-eps_small) = 1;
    f = VOFWLIC_scalarBC2d(f,nx,ny);
    fold = f;    
    
    % Strang splitting, x first if 0, y first if 1
    splitDir = mod(step-1,2);
    % splitDir = 0;
    % remember to keep direction split in column vectors
    if (splitDir == 0) % X -> Y
        % X sweep first
        [normx,normy] = VOFWLIC_normal2d(f,nx,ny,dx,dy);
        for j = 2:ny+1
            fs = f(:,j);
            us = umac(:,j);
            ns = [normx(:,j), normy(:,j)];
            
            fstar = VOFWLIC_flux1d(nx, fs,us,ns, dt,dx, 1,2);
            %
            f(2:nx+1,j) = fstar(2:nx+1) + dt/dx * fs(2:nx+1).*(us(3:nx+2)-us(2:nx+1));
        end
        
        % Y sweep next
        [normx,normy] = VOFWLIC_normal2d(f,nx,ny,dx,dy);
        for i = 2:nx+1
            fs = f(i,:);
            us = vmac(i,:);
            ns = [normx(i,:)', normy(i,:)'];
            
            fstar = VOFWLIC_flux1d(ny, fs',us',ns, dt,dy, 2,2);
            %
            f(i,2:ny+1) = fstar(2:ny+1)' + dt/dy * fs(2:ny+1).*(us(3:ny+2)-us(2:ny+1));
        end
        
    elseif (splitDir == 1) % Y -> X
        % Y sweep first
        [normx,normy] = VOFWLIC_normal2d(f,nx,ny,dx,dy);
        for i = 2:nx+1
            fs = f(i,:);
            us = vmac(i,:);
            ns = [normx(i,:)', normy(i,:)'];
            
            fstar = VOFWLIC_flux1d(ny, fs',us',ns, dt,dy, 2,2);
            %
            f(i,2:ny+1) = fstar(2:ny+1)' + dt/dy * fs(2:ny+1).*(us(3:ny+2)-us(2:ny+1));
        end
        
        % X sweep first
        [normx,normy] = VOFWLIC_normal2d(f,nx,ny,dx,dy);
        for j = 2:ny+1
            fs = f(:,j);
            us = umac(:,j);
            ns = [normx(:,j), normy(:,j)];
            
            fstar = VOFWLIC_flux1d(nx, fs,us,ns, dt,dx, 1,2);
            %
            f(2:nx+1,j) = fstar(2:nx+1) + dt/dx * fs(2:nx+1).*(us(3:nx+2)-us(2:nx+1));
        end
    end
    
    
    if (mod(step,5) == 0)
        fsum = sum(sum(f(2:nx+1,2:ny+1)));
        
        prompt = ['step=',int2str(step),';time=',num2str(time), ...
        ';mass/mass0=',num2str(fsum),'/',num2str(f0sum),'=',num2str(fsum/f0sum)];
        disp(prompt);
        
        subplot(1,2,1);
        contourf(X',Y',f',4); 
        hold on;
        contour(X',Y',f0',[0.5]);
        hold off;
        axis equal;
        axis([0 Lx 0 Ly]);
        title(prompt);
        
        subplot(1,2,2);
        mesh(X',Y',f');
        % surf(X',Y',f','EdgeColor','none');
        axis([0 Lx 0 Ly -0.1 1.1]);
        xlabel('X'); ylabel('Y');
        title(prompt);
        
        drawnow;
    end
end




