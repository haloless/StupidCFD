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

% ## LSTest

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-09

clc;
clear all;

xa = -1; xb = 1;
xlen = xb - xa;
ya = -1; yb = 1;
ylen = yb - ya;

nx = 32;
ny = 32;
hx = xlen / nx;
hy = ylen / ny;

cellxs = linspace(xa-hx/2,xb+hx/2,nx+2);
cellys = linspace(ya-hy/2,yb+hy/2,ny+2);
[X,Y] = ndgrid(cellxs,cellys);

% circle
LSDistFunc = @(x,y) sqrt((x).^2+(y).^2) - 0.25;
% bounding box
LSDistFunc = @(x,y) min(min(abs(x-xa),abs(x-xb)), min(abs(y-ya),abs(y-yb)));
% kink
% LSDistFunc = 
% distorted circle
LSDistFunc = @(x,y) (0.05+(x+1).^2+(y+0.5).^2) .* (sqrt(x.^2+4*y.^2)-1);


Is = 2:nx+1;
Js = 2:ny+1;

ls = LSDistFunc(X,Y);
lsgx = zeros(nx+2,ny+2);
lsgy = zeros(nx+2,ny+2);
lsgx(Is,Js) = 1/(2*hx) * (ls(Is+1,Js)-ls(Is-1,Js));
lsgy(Is,Js) = 1/(2*hy) * (ls(Is,Js+1)-ls(Is,Js-1));

quantity = abs(1 - sqrt(lsgx.^2 + lsgy.^2));
if (1)
    figure;
    [cont,h] = contour(X(Is,Js)',Y(Is,Js)',ls(Is,Js)');
    clabel(cont,h);
    hold on;
    quiver(X(Is,Js)',Y(Is,Js)', lsgx(Is,Js)',lsgy(Is,Js)');
    hold off;
    axis([xa xb ya yb]);
    axis equal;
    
    
    
    figure;
    surf(X(Is,Js)',Y(Is,Js)',quantity(Is,Js)');
    colorbar;
end







