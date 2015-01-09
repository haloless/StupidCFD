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
% ## Created: 2013-07-18

clc;
clear all;

% Description

% constants
ImUnit = 1i; % we keep this for clearity

nu = 1e-3;
% nu = 2e-6;

N = 128; % keep this power of 2
Nx = N;
Ny = N;

Lx = 2*pi;
Ly = 2*pi;
dx = Lx / Nx;
dy = Lx / Ny;

max_time = 100;
max_step = 50000;
dt_max = 1e-2;
cfl = 0.4;
dt = 1e-2;


% initial vorticity
[is,js] = meshgrid(1:Nx,1:Ny);
% [xs,ys] = meshgrid((1:Nx)*dx-dx/2,(1:Ny)*dy-dy/2);
[xs,ys] = meshgrid((1:Nx)*dx-dx,(1:Ny)*dy-dy);

prob_name = 'vortice2';

switch prob_name
    case {'vortice2'}
        sigma = pi/10;
        x1 = 4*pi/5; y1 = x1;
        x2 = 6*pi/5; y2 = x2;
        fw = @(x,y) exp((cos(x-x1)+cos(y-y1)-2) / sigma^2) + exp((cos(x-x2)+cos(y-y2)-2) / sigma^2);
        w0 = fw(xs,ys) - dblquad(fw,0,Lx,0,Ly);
    case {'TaylorGreen'}
        % analytical solution: 
        % w = 2*sin(x)*sin(y)*F(t), F(t) = exp(-2*nu*t)
        w0 = 2 * sin(xs) .* sin(ys);
    otherwise
        error('Unknown prob_name=%s', prob_name);
end

% wave numbers
% XXX: kx is compatible, why ky is different???
% kx = [0, 1, 2, ..., n/2-1, -n/2, -n/2+1, ..., -1]
% ky = [0, -1, -2, ..., -n/2+1, n/2, n/2-1, ..., 1]
kx = ImUnit * ones(1,Ny)' * [0:Nx/2-1,-Nx/2:-1];
ky = ImUnit * [0:-1:-Ny/2+1,Ny/2:-1:1]' * ones(1,Nx);
% 3/2 rule of anti-aliasing
cutoff_filter = abs(kx)<1/3*Nx & abs(ky)<1/3*Ny;
%
k2_visc = kx.^2 + ky.^2;
%
k2_ppe = kx.^2 + ky.^2;
k2_ppe(1,1) = 1; % singular at this point, i.e. kx^2+ky^2 =0
%
w = w0;
w_hat = fft2(w);

istep = 0;
time = 0;

while (time<max_time && istep<max_step)
    
    % % stream function
    % psi_hat = -w_hat ./ k2_ppe;
    
    % % physical velocity
    % u = real(ifft2(ky .* psi_hat));
    % v = real(ifft2(-kx .* psi_hat));
    % dwdx = real(ifft2(kx .* w_hat));
    % dwdy = real(ifft2(ky .* w_hat));
    
    % H = u.*dwdx + v.*dwdy;
    % H_hat = fft2(H);
    % H_hat = cutoff_filter .* H_hat;
    
    % % Crank-Nicholson update
    % visc_coef = 0.5*dt*nu * k2_visc;
    % % visc_coef = visc_coef .* cutoff_filter;
    % w_hat_new = (dt*H_hat + visc_coef.*w_hat + w_hat) ./ (1 - visc_coef);
    
    % TVD-RK3
    L0 = Vorticity_incr(w_hat,kx,ky,nu,k2_visc,k2_ppe,cutoff_filter);
    w_hat1 = w_hat + dt*L0;
    
    L1 = Vorticity_incr(w_hat1,kx,ky,nu,k2_visc,k2_ppe,cutoff_filter);
    w_hat2 = 3/4*w_hat + 1/4*w_hat1 + 1/4*dt*L1;
    
    L2 = Vorticity_incr(w_hat2,kx,ky,nu,k2_visc,k2_ppe,cutoff_filter);
    w_hat_new = 1/3*w_hat + 2/3*w_hat2 + 2/3*dt*L2;
    
    
    time = time + dt;
    istep = istep + 1;
    
    if (mod(istep,100)==0)
        disp(['step=', int2str(istep), ...
            ';time=', num2str(time)]);
    end
    if (mod(istep,500)==0)
        w = real(ifft2(w_hat_new));
        contourf(w,32);
        caxis([min(w0(:)), max(w0(:))]);
        colorbar; shading flat; colormap('jet'); 
        axis equal; 
        % axis([0 Lx 0 Ly]);
        title(num2str(time));
        drawnow;
    end
    
    w_hat = w_hat_new;
end

% post-processing

% draw vorticity
figure;
w = real(ifft2(w_hat));
contourf(w,32);
colorbar; shading flat; colormap('jet');
axis equal; 
% axis([0 Lx 0 Ly]);
title(num2str(time));



