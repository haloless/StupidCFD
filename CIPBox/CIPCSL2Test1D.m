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

% ## CIPCSL2Test1D

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-11-29


clc;
clear all;

xlen = 100.0;
% ncell = 100;
ncell = 300;
hx = xlen / ncell;

ng = 1;
nx = ncell + 2*ng;
lo = 1 + ng;
hi = ncell + ng;
vrange = lo:hi;

cellxs = linspace(-hx/2-hx*(ng-1), xlen+hx/2+hx*(ng-1), nx);
edgexs = linspace(-hx*ng, xlen+hx*ng, nx+1);

% VIA and SIA
fv = zeros(1,nx);
fs = zeros(1,nx+1);

mask = (40<=edgexs & edgexs<=60);
fs(mask) = 1.0;
% fs(mask) = 0.5*sin(2*pi/40 * edgexs(mask));
fv(lo:hi) = 0.5 * (fs(lo:hi) + fs((lo:hi)+1));

% velocity
us = zeros(1,nx+1);
us = 1 + 0.5*sin(2*pi/100 * edgexs);
% us(:) = -1;

if (0)
    plot(edgexs,fs, cellxs,fv, edgexs, us);
    legend('SIA','VIA', 'u-edge');
    return
end

max_time = 100.0;
% max_step = 100000;
max_step = 10 * ncell;
dt = max_time / max_step;
time = 0.0;
step = 0;

while (time<max_time && step<max_step)
    time = time + dt;
    step = step + 1;
    
    if (1)
        % fv(lo-1) = fv(hi);
        % fv(hi+1) = fv(lo);
        for i = 1:ng
            fv(lo-i) = fv(hi-i+1);
            fv(hi+i) = fv(lo+i-1);
        end
        % fs(lo) = fs(hi+1);
        % fs(lo-1) = fs(hi);
        % fs(hi+2) = fs(lo+1);
        for i = 0:ng
            fs(lo-i) = fs(hi+1-i);
            fs(hi+1+i) = fs(lo+i);
        end
    end
    
    fvn = fv;
    fsn = fs;
    flx = fs;
    
    
    cipa = fv;
    cipb = fv;
    cipc = fv;
    % construct CIP
    % NOTE that this must involve ghost cells as well
    irange = lo-1:hi+1;
    cipa(irange) = 1/(hx^2) * (-6*fv(irange) + 3*fs(irange) + 3*fs(irange+1));
    cipb(irange) = 1/(hx) * (6*fv(irange) - 4*fs(irange) - 2*fs(irange+1));
    cipc(irange) = fs(irange);
    
    for i = lo:hi+1
        xx = -us(i) * dt;
        if (us(i) > 0) 
            usgn = 1;
            isup = i-1;
            ivup = i-1;
            xup = hx + xx;
        else
            usgn = -1;
            isup = i+1;
            ivup = i;
            xup = xx;
        end
        
        if(abs(xx)>hx) 
            error('DT too large'); 
        end
        
        a1 = cipa(ivup);
        a2 = cipb(ivup);
        a3 = cipc(ivup);
        
        fstar = a1*(xup^2) + a2*xup + a3;
        
        udiff = usgn * (us(i)-us(isup)) / hx;
        
        % fsn(i) = fstar;
        fsn(i) = fstar - dt*0.5*(fstar+fs(i))*udiff;
        
        %
        flx(i) = a1/3*(xup^3) + a2/2*(xup^2) + a3*xup;
        if (us(i) > 0) 
            flx(i) = fv(ivup)*hx - flx(i);
        else
            flx(i) = -flx(i);
        end
    end
    
    
    % for i = lo:hi
        % fvn(i) = fv(i) + 1/hx*(flx(i)-flx(i+1));
    % end
    fvn(vrange) = fv(vrange) + 1/hx*(flx(vrange)-flx(vrange+1));
    
    fs = fsn;
    fv = fvn;
    
    if (mod(step,10) == 0 || 0)
        mass = sum(fv(vrange))*hx;
        prompt = ['step=',int2str(step), ';time=', num2str(time), ...
            ';mass=',num2str(mass)];
        disp(prompt);
        
        plot(edgexs,fs, cellxs,fv);
        legend('SIA','VIA');
        title(prompt);
        drawnow;
    end
end







