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

% ## SWE1D_Driver

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-02-02


clc;
clear all;


SWE1D_Globals;
GRAV = 9.81;
DENS = 1.0;
HEIGHT_SMALL = 1e-10;
UH1 = 1; UHU1 = 2; UH2 = 3; UHU2 = 4; UBOT = 5;
NUCONS = 5;
QH1 = 1; QVX1 = 2; QH2 = 3; QVX2 = 4; QBOT = 5; 
QGH1 = 6; QC1 = 7; QGH2 = 8; QC2 = 9;
NQPRIM = 9;

BC_PER = 0; BC_NEU = 1; BC_DIR = 2;
bctype = zeros(2,NUCONS); % LOW/HIGH, component
bcfill = zeros(2,NUCONS);
bctype(:,:) = BC_NEU;


NRECONS = 1; % piecewise constant
% NRECONS = 2; % linear reconstruction

% prob = 'InternalDamBreakFlat';
% prob = 'InternalDamBreakFlat2';
prob = 'InternalDamBreakBump';

% ncell = 100;
ncell = 128;
% ncell = 256;
% ncell = 1024;
ng = 4;
nx = ncell + ng*2;
lo = ng + 1;
hi = ng + ncell;
Is = lo:hi;

switch prob
case {'InternalDamBreakFlat'}
    xlo = -5.0; xhi = 5.0; xcen = 0.0;
    maxtime = 1.25;
    DENS = 0.8;
case {'InternalDamBreakFlat2'}
    xlo = -5.0; xhi = 5.0; xcen = 0.0;
    maxtime = 25;
    DENS = 0.99805;
case {'InternalDamBreakBump'}
    xlo = -5.0; xhi = 5.0; xcen = 0.0;
    maxtime = 20;
    DENS = 0.98;
otherwise
    error('Unknown problem')
end
xlen = xhi - xlo;
hx = xlen / ncell;

cellxs = linspace(xlo-hx/2-hx*(ng-1), xhi+hx/2+hx*(ng-1), nx);

% cell state
ucons = zeros(NUCONS,nx);

switch prob
case {'InternalDamBreakFlat'}
    for i = lo:hi
        if (cellxs(i) < xcen)
            ucons(UH1,i) = 0.4;
            ucons(UH2,i) = 0.6;
        else
            ucons(UH1,i) = 0.6;
            ucons(UH2,i) = 0.4;
        end
        ucons(UHU1,i) = 0;
        ucons(UHU2,i) = 0;
        ucons(UBOT,i) = 0;
    end
case {'InternalDamBreakFlat2'}
    for i = lo:hi
        if (cellxs(i) < xcen)
            ucons(UH1,i) = 1.8;
            ucons(UH2,i) = 0.2;
        else
            ucons(UH1,i) = 0.2;
            ucons(UH2,i) = 1.8;
        end
        ucons(UHU1,i) = 0;
        ucons(UHU2,i) = 0;
        ucons(UBOT,i) = -2.0;
    end
case {'InternalDamBreakBump'}
    delta = 0.03;
    eta = 1.0;
    for i = lo:hi
        ucons(UBOT,i) = (0.5-delta) * exp(-cellxs(i)^2);
        if (cellxs(i) < xcen)
            ucons(UH1,i) = 0.5;
        else
            ucons(UH1,i) = delta;
        end
        ucons(UH2,i) = eta - ucons(UH1,i) - ucons(UBOT,i);
        ucons(UHU1,i) = 0;
        ucons(UHU2,i) = 0;
    end
otherwise
    error('Unknown problem')
end

if (DENS > 1)
    error('density ratio > 1');
end

if (0)
    xs = cellxs(Is);
    hs = ucons(UH,Is);
    us = ucons(UHU,Is) ./ hs;
    bs = ucons(UBOT,Is);
    % plot(xs,hs+bs,'x-',xs,bs,'.-', xs,us,'-');
    % legend('surface','bottom', 'velocity');
    plot(xs,hs+bs,'x-',xs,bs,'.-');
    legend('surface','bottom');
    drawnow
end

maxstep = ncell * 100;
% cfl = 0.5;
cfl = 0.8;
dt = maxtime / maxstep;
time = 0.0;
step = 0;

while (time<maxtime && step<maxstep)
    [Ln,dtau] = TwoLayerSWE1D_LuOp(ucons, lo,hi,nx,ng,hx);
    if (1)
        dt = cfl * dtau;
        dt = min([dt, maxtime-time]);
    end
    u1 = ucons + dt*Ln;
    
    [L1,dtau] = TwoLayerSWE1D_LuOp(u1, lo,hi,nx,ng,hx);
    un = 1/2*ucons + 1/2*u1 + 1/2*dt*L1;

    % u2 = 3/4*ucons + 1/4*u1 + 1/4*dt*L1;
    
    % [L2,dtau] = SWE1D_LuOp(u2, lo,hi,nx,ng,hx);
    % un = 1/3*ucons + 2/3*u2 + 2/3*dt*L2;
    
    ucons = un;
    
    time = time + dt;
    step = step + 1;
    
    if (mod(step,10)==0 || time>=maxtime)
        prompt = ['step=',int2str(step),';time=',num2str(time),...
        ';dt=',num2str(dt)];
        disp(prompt);
        
        xs = cellxs(Is);
        h1s = ucons(UH1,Is);
        q1s = ucons(UHU1,Is);
        u1s = q1s ./ h1s;
        h2s = ucons(UH2,Is);
        q2s = ucons(UHU2,Is);
        u2s = q2s ./ h2s;
        bs = ucons(UBOT,Is);
        
        % subplot(2,1,1);
        plot(xs,h1s+h2s+bs,'x-', xs,h2s+bs,'+-', xs,bs,'.-');
        legend('surface 1','surface 2','topography');
        title(prompt);
        
        % subplot(2,1,2);
        % [haxes,hline1,hline2] = plotyy(xs,qs, xs,us);
        % set(hline1,'Marker','.');
        % set(hline2,'Marker','.');
        % legend('discharge','velocity');
        
        drawnow;
    end
end














