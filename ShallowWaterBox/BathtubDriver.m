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

% ## BathtubDriver

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-31

clc;
clear all;

g = 9.8;
u0 = 0;
v0 = 0;
b0 = 0;
h0 = 5000;
h1 = 5030;

% number of nodes
nx = 128 + 1;
ny = 128 + 1;
xa = 0; xb = 100000; xlen = xb - xa;
ya = 0; yb = 100000; ylen = yb - ya;
dx = (xb-xa) / (nx-1);
dy = (yb-ya) / (ny-1);
xs = linspace(xa,xb,nx);
ys = linspace(ya,yb,ny);
[x,y] = ndgrid(xs,ys);

% wave speed
c = sqrt(u0^2+v0^2) + sqrt(g*(h1-b0));

cfl = 0.4;
dt = cfl * min([dx,dy]) / c;

u = zeros(nx,ny);
v = zeros(nx,ny);
b = zeros(nx,ny);

% height
h = zeros(nx,ny);
h(:,:) = h0;
h(find(x>=xa+0.45*xlen & x<=xa+0.55*xlen & y>=ya+0.45*ylen & y<=ya+0.55*ylen)) = h1;
% sea floor
b = zeros(nx,ny);
b(:,:) = b0;
b(find(x<=xa+0.2*xlen)) = h0 / (0.2*xlen) * (0.2*xlen-x(find(x<=xa+0.2*xlen)));

if (0)
    figure; mesh(x,y,h); title('water height');
    figure; mesh(x,y,b); title('floor depth');
end

max_time = 1000;
max_step = 100000;
time = 0.0;
step = 0;

% valid range
I = 2:(nx-1);
J = 2:(ny-1);

while (time<max_time && step<max_step)
    eta = h - b;
    
    hnew = h;
    unew = u;
    vnew = v;
    
    hnew(I,J) = 1/4 * (h(I+1,J) + h(I-1,J) + h(I,J+1) + h(I,J-1)) ...
        - 0.5*dt/dx * u(I,J) .* (eta(I+1,J)-eta(I-1,J)) ...
        - 0.5*dt/dy * v(I,J) .* (eta(I,J+1)-eta(I,J-1)) ...
        - 0.5*dt/dx * (u(I+1,J)-u(I-1,J)) .* eta(I,J) ...
        - 0.5*dt/dy * (v(I,J+1)-v(I,J-1)) .* eta(I,J);
    
    unew(I,J) = 1/4 * (u(I+1,J) + u(I-1,J) + u(I,J+1) + u(I,J-1)) ...
        - 0.5*dt/dx * 1/2 * (u(I+1,J).^2 - u(I-1,J).^2) ...
        - 0.5*dt/dy * v(I,J) * (u(I,J+1)-u(I,J-1)) ...
        - 0.5*dt/dx * g * (h(I+1,J) - h(I-1,J));
    
    vnew(I,J) = 1/4 * (v(I+1,J) + v(I-1,J) + v(I,J+1) + v(I,J-1)) ...
        - 0.5*dt/dy * 1/2 * (v(I,J+1).^2 - v(I,J-1).^2) ...
        - 0.5*dt/dx * u(I,J) * (v(I+1,J)-v(I-1,J)) ...
        - 0.5*dt/dy * g * (h(I,J+1) - h(I,J-1));
    
    % apply BC
    unew(1,:) = 2.5*unew(2,:) - 2*unew(3,:) + 0.5*unew(4,:);
    unew(nx,:) = 2.5*unew(nx-1,:) - 2*unew(nx-2,:) + 0.5*unew(nx-3,:);
    unew(:,1) = 2.5*unew(:,2) - 2*unew(:,3) + 0.5*unew(:,4);
    unew(:,ny) = 2.5*unew(:,ny-1) - 2*unew(:,ny-2) + 0.5*unew(:,ny-3);
    %
    vnew(1,:) = 2.5*vnew(2,:) - 2*vnew(3,:) + 0.5*vnew(4,:);
    vnew(nx,:) = 2.5*vnew(nx-1,:) - 2*vnew(nx-2,:) + 0.5*vnew(nx-3,:);
    vnew(:,1) = 2.5*vnew(:,2) - 2*vnew(:,3) + 0.5*vnew(:,4);
    vnew(:,ny) = 2.5*vnew(:,ny-1) - 2*vnew(:,ny-2) + 0.5*vnew(:,ny-3);
    %
    hnew(1,:) = 2.5*hnew(2,:) - 2*hnew(3,:) + 0.5*hnew(4,:);
    hnew(nx,:) = 2.5*hnew(nx-1,:) - 2*hnew(nx-2,:) + 0.5*hnew(nx-3,:);
    hnew(:,1) = 2.5*hnew(:,2) - 2*hnew(:,3) + 0.5*hnew(:,4);
    hnew(:,ny) = 2.5*hnew(:,ny-1) - 2*hnew(:,ny-2) + 0.5*hnew(:,ny-3);
    
    h = hnew;
    u = unew;
    v = vnew;
    
    
    step = step + 1;
    time = time + dt;
    
    if (mod(step,10)==0)
        disp(['step=',int2str(step),';time=',num2str(time)]);
        mesh(x,y,h);
        % hold on;
        % mesh(x,y,b);
        % hold off;
        axis([xa xb ya yb h0*0.999 h0*1.001]);
        xlabel('x'); ylabel('y');
        zlabel('height');
        drawnow;
    end
    
end




