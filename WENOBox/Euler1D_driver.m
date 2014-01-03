% ## Copyright (C) 2014 homu
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

% ## Euler1D_driver

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-01-03

clc;
clear all;

% global values
Euler1D_globals;
% euler_gamma = 1.4;
GAMMA = 1.4;
URHO = 1; UMX = 2; UETOT = 3;
QRHO = 1; QVX = 2; QPRES = 3;

prob = 'SodShocktube';
% prob = 'LaxTest';

switch prob
case {'SodShocktube'}
    xlo = 0.0;
    xhi = 1.0;
case {'LaxTest'}
    xlo = -5.0;
    xhi = 5.0;
otherwise
    error('Unknown problem');
end
xlen = xhi - xlo;
ncell = 128;
hx = xlen / ncell;

ng = 4;
nx = ncell + 2*ng;
lo = 1 + ng;
hi = ncell + ng;
vrange = lo:hi;

cellxs = linspace(xlo-hx/2-hx*(ng-1), xhi+hx/2+hx*(ng-1), nx);
edgexs = linspace(xlo-hx*ng, xhi+hx*ng, nx+1);

% cell state and flux
ucons = zeros(3,nx);
uprim = zeros(3,nx);

switch prob
case {'SodShocktube'} % Sod's shock tube
    for i = lo:hi
        if(cellxs(i) < 0.5)
            uprim(:,i) = [1.0; 0.0; 1.0];
        else
            uprim(:,i) = [0.125; 0.0; 0.1];
        end
    end
    max_time = 0.2;
case {'LaxTest'} % Lax's Riemann problem
    for i = lo:hi
        if(cellxs(i) < 0.0)
            uprim(:,i) = [0.445; 0.698; 3.528];
        else
            uprim(:,i) = [0.5; 0.0; 0.571];
        end
    end
    max_time = 1.3;
otherwise
    error('Unknown problem.')
end
ucons = Euler1D_PrimToCons(uprim,ucons,lo,hi);


max_step = 100 * ncell;
cfl = 0.5;
dt = max_time / max_step;
time = 0.0;
step = 0;

while (time<max_time && step<max_step)
    
    [Ln,dtau] = Euler1D_LuOp(ucons, lo,hi,nx,ng,hx);
    if (1)
        dt = cfl * dtau;
        dt = min([dt, max_time-time]);
    end
    u1 = ucons + dt*Ln;
    
    [L1,dtau] = Euler1D_LuOp(u1, lo,hi,nx,ng,hx);
    u2 = 3/4*ucons + 1/4*u1 + 1/4*dt*L1;
    
    [L2,dtau] = Euler1D_LuOp(u2, lo,hi,nx,ng,hx);
    un = 1/3*ucons + 2/3*u2 + 2/3*dt*L2;
    
    ucons = un;
    
    
    time = time + dt;
    step = step + 1;
    
    
    if (mod(step,10)==0 || time>=max_time)
        prompt = ['step=',int2str(step),';time=',num2str(time),';dt=',num2str(dt)];
        disp(prompt);
        
        % conservative -> primitive
        uprim = Euler1D_ConsToPrim(ucons,uprim, lo,hi);
        xs = cellxs(vrange);
        plot(xs,uprim(QRHO,vrange),'x-', xs,uprim(QVX,vrange),'x-', xs,uprim(QPRES,vrange),'x-');
        legend('density','velocity','pressure');
        title(prompt);
        drawnow;
    end
end
















