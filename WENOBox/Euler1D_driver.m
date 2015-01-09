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
GAMMA = 1.4;
URHO = 1; UMX = 2; UETOT = 3;
QRHO = 1; QVX = 2; QPRES = 3;

% boundary condition
BC_TRANSPARENT = 0;
BC_REFLECT = 1;
bctype = zeros(2); bctype(:) = BC_TRANSPARENT;

RHO_SMALL = 1e-13;

% do_positivity_limit = 1;
do_positivity_limit = 0;

% prob = 'SodShocktube';
% prob = 'LaxTest';
prob = 'DoubleRarefaction';
% prob = 'LeftBlastWave';
% prob = 'SedovBlastWave';
% prob = 'WoodwardCollelaBlastWave';
% prob = 'ShuOsherShockEntropy';
% prob = 'LeblancShocktube';



switch prob
case {'SodShocktube'}
    xlo = 0.0; xhi = 1.0;
case {'LaxTest'}
    xlo = -5.0; xhi = 5.0;
case {'DoubleRarefaction'}
    % xlo = -1.0; xhi = 1.0;
    xlo = -0.5; xhi = 0.5;
case {'LeftBlastWave'}
    xlo = 0.0; xhi = 1.0;
case {'SedovBlastWave'}
    xlo = -2.0; xhi = 2.0;
case {'WoodwardCollelaBlastWave'}
    xlo = 0.0; xhi = 1.0;
    bctype(:) = BC_REFLECT;
case {'ShuOsherShockEntropy'}
    xlo = -5.0; xhi = 5.0;
case {'LeblancShocktube'}
    xlo = -10.0; xhi = 10.0;
otherwise
    error('Unknown problem');
end
xlen = xhi - xlo;
ncell = 128;
% ncell = 256;
% ncell = 512;
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
case {'DoubleRarefaction'}
    for i = lo:hi
        if (cellxs(i) < 0)
            % uprim(:,i) = [7.0; -1.0; 0.2];
            uprim(:,i) = [1.0; -2.0; 0.4];
        else
            % uprim(:,i) = [7.0; 1.0; 0.2];
            uprim(:,i) = [1.0; 2.0; 0.4];
        end
    end
    % max_time = 0.6;
    max_time = 0.15;
case {'LeftBlastWave'}
    for i = lo:hi
        if (cellxs(i) < 0.5)
            uprim(:,i) = [1.0; 0.0; 1000.0];
        else
            uprim(:,i) = [1.0; 0.0; 0.01];
        end
    end
    max_time = 0.012;
case {'SedovBlastWave'}
    for i = lo:hi
        if (abs(cellxs(i)) < hx)
            uprim(:,i) = [1.0; 0.0; (GAMMA-1)*3200000/hx];
        else
            uprim(:,i) = [1.0; 0.0; 1e-12];
        end
    end
    max_time = 0.001;
case {'WoodwardCollelaBlastWave'}
    for i = lo:hi
        if (cellxs(i) < 0.1) % left
            uprim(:,i) = [1.0; 0.0; 1000.0];
        elseif (cellxs(i) >= 0.9) % right
            uprim(:,i) = [1.0; 0.0; 100.0];
        else % middle
            uprim(:,i) = [1.0; 0.0; 0.01];
        end
    end
    max_time = 0.038;
case {'ShuOsherShockEntropy'}
    for i = lo:hi
        if (cellxs(i) < -4.0)
            uprim(:,i) = [3.857143; 2.629369; 10.33333];
        else
            uprim(:,i) = [1+0.2*sin(5*cellxs(i)); 0.0; 1.0];
        end
    end
    max_time = 1.8;
case {'LeblancShocktube'}
    for i = lo:hi
        if (cellxs(i) < 0.0)
            uprim(:,i) = [2.0; 0.0; 1.0e9];
        else
            uprim(:,i) = [0.001; 0.0; 1.0];
        end
    end
    max_time = 0.0001;
otherwise
    error('Unknown problem.')
end
ucons = Euler1D_PrimToCons(uprim,ucons,lo,hi);


max_step = 100 * ncell;
% NOTE CFL must be set small if positivity limiter is used.
cfl = 0.5;
% cfl = 0.1;
% cfl = 0.05;
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
        eint = 1/(GAMMA-1) .* uprim(QPRES,vrange) ./ uprim(QRHO,vrange);
        xs = cellxs(vrange);
        plot(xs,uprim(QRHO,vrange),'x-', xs,uprim(QVX,vrange),'x-', ...
        xs,uprim(QPRES,vrange),'x-', xs, eint);
        legend('density','velocity','pressure','internal');
        title(prompt);
        drawnow;
    end
end
















