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

% ## WENO_BurgersDriver

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-12

% Description:
% du/dt + dF(u)/dx = 0

clc;
clear all;

% global variables
WENO_BurgersGlobals;

% flux F = F(u), must support vectorized operations
probno = 2;
switch probno
    case {1}
        % linear scalar transport
        % F = a * u
        a = -0.5;
        F = @(u) a * u;
        dFdu = @(u) a * ones(size(u));
    case {2}
        % Burgers equation
        % F = 1/2 * u^2
        F = @(u) 0.5 * u.^2;
        dFdu = @(u) u;
    case {3}
        % nonlinear nonconvex scalar Buckley-Leverett
        F = @(u) (4*u.^2) ./ (4*u.^2 + (1-u).^2);
        dFdu = @(u) (8*u.*(1-u)) ./ (4*u.^2 + (1-u).^2).^2;
    otherwise
        error('Unknown prob_no=%d', probno);
end

cfl = 0.25;

x_lo = 0;
x_hi = 2;
ncell = 80;
nnode = ncell + 2;
dx = (x_hi-x_lo) / ncell;

% WENO
r = 3;
k = 2;
ngrow = k;
nx = nnode + ngrow*2;
xs = linspace(x_lo-dx/2-ngrow*dx,x_hi+dx/2+ngrow*dx,nx);
valid_range = 1+ngrow:nx-ngrow;

% initial state
u0 = 0.5 + sin(1*pi*xs);
% u0 = xs>=0.5 & xs<=1;

u = u0;
u_new = zeros(size(u));

time = 0;
% max_time = 1.0;
max_time = 1.5/pi;
istep = 0;
max_step = 5000;

figure;

while (time<max_time & istep<=max_step)
    dt = cfl * dx / max(abs(u));
    dt = min([dt,max_time-time]);
    
    % % apply BC
    % u = WENO_bc(u, nx, ngrow);
    
    % % flux splitting
    % [fp,fn] = WENO_fluxsplit(u,F,dFdu);
    
    % % reconstruct flux at cells interface
    % hp = zeros(size(fp));
    % hn = zeros(size(fn));
    % for i = valid_range
        % xr = i-ngrow:i+ngrow;
        % [h_right,h_left] = WENO_flux(fp(xr),fn(xr));
        % hn(i) = h_right;
        % hp(i-1) = h_left;
    % end
    
    % % f_1/2 = f_1/2^plus + f_1/2^minus
    % h = hp + hn;
    
    % % h_right = h;
    % % h_left = [0, h(1:end-1)];
    
    % % update
    % % u_new = u - dt/dx * (h_right-h_left);
    
    % u_new(2:end) = u(2:end) - dt/dx * diff(h);
    
    % L = WENO_LuOp(u, F, dFdu, nx, ngrow, dx);
    % u_new = u + dt*L;
    
    % TVD-RK3 time stepping
    Ln = WENO_LuOp(u, F, dFdu, nx, ngrow, dx);
    u1 = u + dt*Ln;
    
    L1 = WENO_LuOp(u1, F, dFdu, nx, ngrow, dx);
    u2 = 3/4*u + 1/4*u1 + 1/4*dt*L1;
    
    L2 = WENO_LuOp(u2, F, dFdu, nx, ngrow, dx);
    u_new = 1/3*u + 2/3*u2 + 2/3*dt*L2;
    
    u_new = WENO_bc(u_new, nx, ngrow);
    u = u_new;
    
    time = time + dt;    
    istep = istep + 1;
    if(mod(istep,10)==0)
        disp(['step=',int2str(istep),';time=',num2str(time),';dt=',num2str(dt)]);
        
        plot(xs(valid_range),u(valid_range),'-');
        axis([x_lo x_hi min(u0)-0.1 max(u0)+0.1]);
        % grid on;
        drawnow;
    end
    
end



