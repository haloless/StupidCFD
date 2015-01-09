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

% ## RiderKotheVortexTracing

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-10-16

clc;
clear all;

xlen = 1.0;
ylen = 1.0;

% Rider-Kothe single vortex parameters
cx = xlen * 0.5;
cy = ylen * 0.75;
cr = xlen * 0.15;
T = 8.0;

% generate points for tracing
% np = 512;
np = 1024;
theta = (2*pi/np) * (1:np)';
xs0 = cx + cr*cos(theta);
ys0 = cy + cr*sin(theta);

if (1)
    figure;
    plot(xs0,ys0,'-');
    axis equal;
    axis([0 xlen 0 ylen]);
end



max_step = 4000;
max_time = T * 0.5;
%max_step = 8000
%max_time = T
dtau = max_time / max_step;
time = 0;

xs = xs0;
ys = ys0;
[us,vs] = RiderKotheVortex(xs,ys,time,T);

% TODO change this to a method-of-line
for step = 1:max_step
    
    
    % update position
    dtau_half = dtau * 0.5;
    % RK4
    [u1,v1] = RiderKotheVortex(xs,ys,time,T);
    [u2,v2] = RiderKotheVortex(xs+dtau_half*u1,ys+dtau_half*v1, time+dtau_half, T);
    [u3,v3] = RiderKotheVortex(xs+dtau_half*u2,ys+dtau_half*v2, time+dtau_half, T);
    [u4,v4] = RiderKotheVortex(xs+dtau*u3,ys+dtau*v3, time+dtau,T);
    
    xs = xs + dtau/6 * (u1 + 2*u2 + 2*u3 + u4);
    ys = ys + dtau/6 * (v1 + 2*v2 + 2*v3 + v4);
    
    time = time + dtau;
    
    if (mod(step,50) == 0) 
        prompt = ['step=',int2str(step),';time=',num2str(time)];
        disp(prompt);
        
        plot(xs0,ys0,'k--', [xs;xs(1)],[ys;ys(1)],'r-');
        axis equal;
        axis([0 xlen 0 ylen]);
        %title(prompt);
        drawnow;
    end
end




