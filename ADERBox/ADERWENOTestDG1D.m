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

% ## ADERWENOTestDG1D

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-01-20



% clc
clear all

% ader_order = 3;
ader_order = 4;
% ader_order = 5;

ADERWENOGlobals1D;
ADERWENOInit1D(ader_order);

M = MDegree;
N = NPoint;
Nd = N * N;

eta = GausEta;
wgt = GausWgt;

use_modal = 0;
% use_modal = 1;
SpaceTimeDGGlobals1D;
SpaceTimeDGInit1D(1,use_modal);

% test
% solve a simple linear ODE
% du/dt = -nu*u
% nu = 3;
nu = 10;
% nu = 100;
% nu = 1000;
u0 = 1.0;
u = u0;
uimp = u0;
ucn = u0;
% nstep = 100;
nstep = 25;
dt = 1.0 / nstep;
dx = 1.0; % fake DX

us = zeros(M*2+1,1);
results = [0, u, uimp, ucn];

time = 0;
% maxstep = 1;
maxstep = nstep;
for step = 1:maxstep
    time = time + dt;
    
    nustar = nu * dt;
    astar = 0;
    
    % fake the spatial reconstruction
    if (use_modal)
        w0 = zeros(N,1);
        w0(1) = u;
    else
        % w0 = zeros(N,1);
        % w0(:) = u;
        
        us(:) = u;
        [w0, p0] = ADERWENOReconstruct1D(us,dx);
    end
    
    A = F1Mat - KtMat + astar*KxMat + nustar*MMat;
    % A = F1Mat - KtMat + nu*MMat;
    % A(abs(A)<1e-10) = 0;
    rhs = F0Mat * w0;
    uh = A \ rhs;
    
    sbar = 0;
    swgt = 0;
    
    for p = 1:Nd
        [px,pt] = ind2sub([N N], p);
        if(use_modal)
            vpx = RescaledLegendrePolyVal(px-1,eta);
            vpt = RescaledLegendrePolyVal(pt-1,eta);
            sbar = sbar + uh(p)*sum(vpx.*wgt)*sum(vpt.*wgt);
        else
            if (0)
                vpx = polyval(LagrPsi(:,px),eta);
                vpt = polyval(LagrPsi(:,pt),eta);
                sbar = sbar + uh(p)*sum(vpx.*wgt)*sum(vpt.*wgt);
                swgt = swgt + sum(vpx.*wgt)*sum(vpt.*wgt);
            else
                sbar = sbar + uh(p)*wgt(px)*wgt(pt);
                swgt = swgt + wgt(px)*wgt(pt);
            end
        end
    end
    
    unew = u - dt*nu*sbar;
    u = unew;
    
    if (1)
        prompt = ['step=',int2str(step),';time=',num2str(time)];
        disp(prompt);
    end
    
    % simple implicit solution
    uimp1 = uimp / (1+dt*nu);
    % Crank-Nicolson
    ucn1 = (1-dt*nu/2)*ucn / (1+dt*nu/2);
    uimp = uimp1;
    ucn = ucn1;
    results(end+1,:) = [time,u,uimp,ucn];
end

if (1)
    ts = results(:,1);
    
    figure;
    plot(ts,u0*exp(-nu*ts),'-', ts,results(:,2),'x', ...
    ts,results(:,3),'o', ts,results(:,4),'s');
    % semilogy(ts,u0*exp(-nu*ts),'-', ts,results(:,2),'.');
    legend('ana','DG','Implicit','Crank-Nicolson');
end




























