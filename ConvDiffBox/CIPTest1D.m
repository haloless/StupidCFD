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

% ## CIPTest1D

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-05

clc;
clear all;

xlen = 1.0;

ncell = 256;
hx = xlen / ncell;

nx = ncell + 2;
xs = linspace(-hx/2,xlen+hx/2,nx);

c = -1.0;
us = c * ones(1,nx);
fs = zeros(1,nx);
gs = zeros(1,nx);

do_tan_tx = 1;
tan_coef = 0.9;
tan_trans = @(f) tan((f-0.5)*tan_coef*pi);
tan_itrans = @(tf) atan(tf) / (tan_coef*pi) + 0.5;

prob = 'SquaredWave';
switch prob
% case {'SmoothGaussian'}
    % f0 = exp(-8 * (2*));
case {'SquaredWave'}
    fs0 = 1.0 * (0.1<=xs & xs<0.2);
    gs0 = zeros(1,nx);
    gs0(2:ncell+1) = 1/(2*hx) * (fs0(3:ncell+2)-fs0(1:ncell));
    if (do_tan_tx)
        fs0t = tan_trans(fs0);
        gs0(2:ncell+1) = 1/(2*hx) * (fs0t(3:ncell+2)-fs0t(1:ncell));
    end
    % gs0(2:ncell+1) = -1/hx * (fs0(1:ncell)-fs0(2:ncell+1));
otherwise
    error('Unknown problem: %s',prob);
end

fs = fs0;
gs = gs0;

if (0)
    figure;
    subplot(2,1,1); plot(xs,fs, xs,gs, xs,us); legend('f0','g0','u');
    subplot(2,1,2); plot(xs,fs, xs,us); legend('f0','u');
end

cfl = 0.2;
dt = cfl * hx / max(abs(us));
max_time = 1.5;
max_step = 10000;
time = 0;
step = 0;

figure;
while (time<max_time && step<max_step)
    time = time + dt;
    step = step + 1;
    
    % BC
    fs(1) = fs(ncell+1);
    fs(ncell+2) = fs(2);
    gs(1) = gs(ncell+1);
    gs(ncell+2) = gs(2);
    
    if (do_tan_tx)
        fs = tan_trans(fs);
    end
    
    
    fn = fs;
    gn = gs;
    
    for i = 2:ncell+1
        if (us(i)>0)
            isign = -1;
        else
            isign = 1;
        end
        iup = i + isign;
        
        xx = -us(i) * dt;
        fdif = isign * (fs(iup)-fs(i)) / hx;
        xam1 = (gs(i)+gs(iup)-2*fdif) / hx^2;
        xbm1 = isign * (3*fdif - 2*gs(i) - gs(iup)) / hx;
        
        fn(i) = xam1*xx^3 + xbm1*xx^2 + gs(i)*xx + fs(i);
        gn(i) = 3*xam1*xx^2 + 2*xbm1*xx + gs(i);
        % fn(i) = ((xam1*xx+xbm1)*xx + gs(i)) * xx + fs(i);
        % gn(i) = (3*xam1*xx + 2*xbm1) * xx + gs(i);
    end
    
    if (do_tan_tx)
        fn = tan_itrans(fn);
    end
    
    fs = fn;
    gs = gn;
    
    if (mod(step,20)==0)
        valid = 2:ncell+1;
        xadv = us(valid) * time;
        xadv = xadv - xlen*floor(xadv/xlen);
        xadv = xs(valid) + xadv;
        fadv = fs0(valid);
        il = find(xadv<0); ir = find(xadv>=xlen);
        xadv(il) = xadv(il) + xlen;
        xadv(ir) = xadv(ir) - xlen;
        [xadv,ix] = sort(xadv);
        fadv = fadv(ix);
        
        plot(xadv,fadv,'-', xs(valid),fs(valid),'o-'); legend('f-exact','f-CIP');
        title([num2str(time)]);
        % pause(1/20);
        drawnow;
    end
    
end







