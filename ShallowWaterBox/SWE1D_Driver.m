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
GRAV = 9.812;
HEIGHT_SMALL = 1e-10;
UH = 1; UHU = 2; UBOT = 3;
NUCONS = 3;
QH = 1; QVX = 2; QBOT = 3; QGH = 4; QC = 5;
NQPRIM = 5;

BC_PER = 0; BC_NEU = 1; BC_DIR = 2;
bctype = zeros(2,NUCONS); % LOW/HIGH, component
bcfill = zeros(2,NUCONS);
bctype(:,:) = BC_NEU;


NRECONS = 1; % piecewise constant
NRECONS = 2; % linear reconstruction
% NRECONS = 3; % ADER-WENO 3
% NRECONS = 5; % ADER-WENO 5

if (NRECONS >= 3)
    ADERWENOGlobals1D;
    ADERWENOInit1D(NRECONS);
end


% prob = 'RP1';
% prob = 'RP2'; % dam break with flat bottom
% prob = 'RP3'; % dam break with wet and step deck
% prob = 'RP4';
% prob = 'RP5';
% prob = 'RP6';
% prob = 'Bump_a'; % transcritical flow without shock
% prob = 'Bump_b'; % transcritical flow with shock
prob = 'Bump_c'; % subcritical


ncell = 100;
% ncell = 200;
% ncell = 128;
% ncell = 256;
% ncell = 1024;
ng = 8;
nx = ncell + ng*2;
lo = ng + 1;
hi = ng + ncell;
Is = lo:hi;

switch prob
case {'RP1'}
    xlo = 0.0; xhi = 1.0; xcen = 0.5;
    maxtime = 0.1;
case {'RP2'}
    xlo = 0.0; xhi = 1.0; xcen = 0.5;
    maxtime = 0.075;
case {'RP3'}
    xlo = -5.0; xhi = 5.0; xcen = 0.0;
    maxtime = 1.0;
case {'RP4'}
    xlo = -3.0; xhi = 4.0; xcen = 0.0;
    maxtime = 1.0;
case {'RP5'}
    xlo = -15.0; xhi = 5.0; xcen = 0.0;
    maxtime = 1.0;
case {'RP6'}
    xlo = -10.0; xhi = 4.0; xcen = 0.0;
    maxtime = 1.0;
case {'Bump_a','Bump_b','Bump_c'}
    xlo = 0.0; xhi = 20.0; xcen = 10.0;
    maxtime = 200;
otherwise
    error('Unknown problem')
end
xlen = xhi - xlo;
hx = xlen / ncell;

cellxs = linspace(xlo-hx/2-hx*(ng-1), xhi+hx/2+hx*(ng-1), nx);

% cell state and flux
ucons = zeros(NUCONS,nx);
uprim = zeros(NQPRIM,nx);

switch prob
case {'RP1'}
    for i = lo:hi
        if (cellxs(i) < xcen)
            uprim(QH,i) = 2.0; uprim(QVX,i) = 0.0; uprim(QBOT,i) = 0.0;
        else
            uprim(QH,i) = 1.0; uprim(QVX,i) = 0.0; uprim(QBOT,i) = 1.0;
        end
    end
case {'RP2'} 
    for i = lo:hi
        if (cellxs(i) < xcen)
            uprim(QH,i) = 1.0; uprim(QVX,i) = 0.0; uprim(QBOT,i) = 0.0;
        else
            uprim(QH,i) = HEIGHT_SMALL; uprim(QVX,i) = 0.0; uprim(QBOT,i) = 0.0;
        end
    end
case {'RP3'}
    for i = lo:hi
        if (cellxs(i) < xcen)
            uprim(QH,i) = 1.46184; uprim(QVX,i) = 0.0; uprim(QBOT,i) = 0.0;
        else
            uprim(QH,i) = 0.30873; uprim(QVX,i) = 0.0; uprim(QBOT,i) = 0.2;
        end
    end
case {'RP4'}
    for i = lo:hi
        if (cellxs(i) < xcen)
            uprim(QH,i) = 0.56900; uprim(QVX,i) = 2.12634; uprim(QBOT,i) = 0.0;
        else
            uprim(QH,i) = 0.60805; uprim(QVX,i) = 0.0; uprim(QBOT,i) = 0.2;
        end
    end
case {'RP5'}
    for i = lo:hi
        if (cellxs(i) < xcen)
            uprim(QH,i) = 0.75; uprim(QVX,i) = -9.49365; uprim(QBOT,i) = 0.0;
        else
            uprim(QH,i) = 1.10594; uprim(QVX,i) = -4.94074; uprim(QBOT,i) = 0.2;
        end
    end
case {'RP6'}
    for i = lo:hi
        if (cellxs(i) < xcen)
            uprim(QH,i) = 0.75; uprim(QVX,i) = -1.35624; uprim(QBOT,i) = 0.0;
        else
            uprim(QH,i) = 1.10594; uprim(QVX,i) = -4.94074; uprim(QBOT,i) = 0.2;
        end
    end
case {'Bump_a','Bump_b','Bump_c'}
    for i = lo:hi
        if (cellxs(i)>=8.0 && cellxs(i)<=12.0)
            uprim(QBOT,i) = 0.2 - 0.05*(cellxs(i)-10)^2;
        else
            uprim(QBOT,i) = 0;
        end
        % uprim(QH,i) = 0.33 - uprim(QBOT,i);
        uprim(QH,i) = 0.5 - uprim(QBOT,i);
        uprim(QVX,i) = 0;
    end
    if (strcmp(prob, 'Bump_a'))
        qin = 1.53; hout = 0.66;
    elseif (strcmp(prob, 'Bump_b'))
        qin = 0.18; hout = 0.33;
    elseif (strcmp(prob, 'Bump_c'))
        qin = 4.42; hout = 2.0;
    end
    
    bctype(1,UH) = BC_NEU; 
    bctype(2,UH) = BC_DIR; bcfill(2,UH) = hout;
    bctype(1,UHU) = BC_DIR; bcfill(1,UHU) = qin;
    bctype(2,UHU) = BC_NEU;
otherwise
    error('Unknown problem')
end

ucons(:,lo:hi) = SWE1D_PrimToCons(uprim(:,lo:hi));

if (1)
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
    [Ln,dtau] = SWE1D_LuOp(ucons, lo,hi,nx,ng,hx);
    if (1)
        dt = cfl * dtau;
        dt = min([dt, maxtime-time]);
    end
    u1 = ucons + dt*Ln;
    un = u1;
    if(1)
    [L1,dtau] = SWE1D_LuOp(u1, lo,hi,nx,ng,hx);
    un = 1/2*ucons + 1/2*u1 + 1/2*dt*L1;
    end

    % u2 = 3/4*ucons + 1/4*u1 + 1/4*dt*L1;
    
    % [L2,dtau] = SWE1D_LuOp(u2, lo,hi,nx,ng,hx);
    % un = 1/3*ucons + 2/3*u2 + 2/3*dt*L2;
    
    ucons = un;
    
    time = time + dt;
    step = step + 1;
    
    if (mod(step,10)==0 || time>=maxtime)
        prompt = ['step=',int2str(step),';time=',num2str(time),...
        ';order=',int2str(NRECONS),';dt=',num2str(dt)];
        disp(prompt);
        
        % uprim = SWE1D_Cons2Prim(ucons);
        xs = cellxs(Is);
        hs = ucons(UH,Is);
        qs = ucons(UHU,Is);
        us = qs ./ hs;
        bs = ucons(UBOT,Is);
        
        % subplot(2,1,1);
        plot(xs,hs+bs,'x-',xs,bs,'.-');
        legend('surface','topography');
        title(prompt);
        
        % subplot(2,1,2);
        % [haxes,hline1,hline2] = plotyy(xs,qs, xs,us);
        % set(hline1,'Marker','.');
        % set(hline2,'Marker','.');
        % legend('discharge','velocity');
        
        drawnow;
    end
end














