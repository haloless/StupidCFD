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

ng = 2;
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
% us(:) = 1;

if (0)
    plot(edgexs,fs, cellxs,fv, edgexs, us);
    legend('SIA','VIA', 'u-edge');
    return
end

max_time = 100.0;
% max_step = 100000;
max_step = 10 * ncell;
dt = max_time / max_step;
% max_step = 1;
% dt = 0.01;
time = 0.0;
step = 0;

while (time<max_time && step<max_step)
    time = time + dt;
    step = step + 1;
    
    if (1)
        for i = 1:ng
            fv(lo-i) = fv(hi-i+1);
            fv(hi+i) = fv(lo+i-1);
        end
        for i = 0:ng
            fs(lo-i) = fs(hi+1-i);
            fs(hi+1+i) = fs(lo+i);
        end
    end
    
    fvn = fv;
    fsn = fs;
    
    fc = fv; % interpolated value at cell center
    dc = fv; % limited slope

    % construct CIP
    % cell center value
    fc(lo-2:hi+2) = 1.5*fv(lo-2:hi+2) - 0.25*fs(lo-2:hi+2) - 0.25*fs(lo-1:hi+3);
    
    % NOTE that this must involve ghost cells as well
    irange = (lo-1:hi+1);
    % limited slope
    slft = fv; slft(irange) = fc(irange) - fc(irange-1);
    srgt = fv; srgt(irange) = fc(irange+1) - fc(irange);
    scen = fv; scen(irange) = 0.5*(fc(irange+1)-fc(irange-1));
    dc(:) = 0.0;
    mask = (slft.*srgt > 0);
    dc(mask) = sign(scen(mask)) .* min(abs(scen(mask)),2*min(abs(slft(mask)),abs(srgt(mask))));
    % for i = lo-1:hi+1
        % if (slft(i)*srgt(i) > 0)
            % smin = min(2*abs(slft(i)), 2*abs(srgt(i)));
            % smin = min(abs(scen(i)), smin);
            % dc(i) = sign(slft(i)) * smin;
        % else
            % dc(i) = 0;
        % end
    % end
    
    c0 = fv;
    c1 = fv;
    c2 = fv;
    c3 = fv;
    
    c0(irange) = fs(irange);
    c1(irange) = 1/(hx) * (6*(fv(irange)-fs(irange)) - 2*dc(irange));
    c2(irange) = 1/(hx^2) * (-3*fs(irange+1) + 9*fs(irange) - 6*fv(irange) + 6*dc(irange));
    c3(irange) = 1/(hx^3) * (4*(fs(irange+1)-fs(irange)) - 4*dc(irange));
    
    % conservative flux
    flx = fs;
    
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
        
        a0 = c0(ivup);
        a1 = c1(ivup);
        a2 = c2(ivup);
        a3 = c3(ivup);
        
        fstar = a0 + a1*(xup) + a2*(xup^2) + a3*(xup^3);
        
        udiff = usgn * (us(i)-us(isup)) / hx;
        
        fsn(i) = fstar - dt*0.5*(fstar+fs(i))*udiff;
        
        %
        flx(i) = a0*xup + a1/2*(xup^2) + a2/3*(xup^3) + a3/4*(xup^4);
        if (us(i) > 0) 
            flx(i) = fv(ivup)*hx - flx(i);
        else
            flx(i) = -flx(i);
        end
    end
    
    % conservative update VIA
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







