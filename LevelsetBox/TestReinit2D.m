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

% ## TestReinit2D

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-04-28

clc;
clear all;

LSGlobals2D;
LSGlobalsInit2D();

%
% dens_spread = 2;
% dhdiv = 1.5;
% dhdiv = 4;

% nx = 32; ny = 32;
nx = 50; ny = 50;
ng = dens_spread + 1;

lo = [ng+1, ng+1];
hi = [nx+ng, ny+ng];

% xlo = -1;
xlo = 0;
xhi = 1;
% ylo = -1;
ylo = 0;
yhi = 1;

dx = (xhi-xlo) / nx;
dy = (yhi-ylo) / ny;
dh = dx;

xpos = linspace(xlo-dx*ng, xhi+dx*ng, nx+ng*2);
ypos = linspace(ylo-dy*ng, yhi+dy*ng, ny+ng*2);
[xpos,ypos] = ndgrid(xpos,ypos);

%% circle
% d0 = xpos.^2 + ypos.^2 - 0.5^2;
% dext = sqrt(xpos.^2 + ypos.^2) - 0.5;

%% line
% d0 = zeros(nx+ng*2,ny+ng*2);
% d0(xpos>0.5) = 0.5 * dh;
% d0(xpos<=0.5) = -0.5 * dh;
% dext = xpos - 0.5;

%% square
square_xlo = 0.0;
square_xhi = 0.3;
square_ylo = 0.0;
square_yhi = 0.4;
d0 = zeros(nx+ng*2,ny+ng*2);
mask = (square_xlo<xpos & xpos<square_xhi & square_ylo<ypos & ypos<square_yhi);
d0(mask) = (2*1.0-1)*0.5*dh;
d0(~mask) = (2*0.0-1)*0.5*dh;
dext = d0;
for j = lo(2):hi(2)
for i = lo(1):hi(1)
    if dext(i,j) > 0
        dist = 99999;
        dist = min(dist, xpos(i,j)-square_xlo);
        dist = min(dist, square_xhi-xpos(i,j));
        dist = min(dist, ypos(i,j)-square_ylo);
        dist = min(dist, square_yhi-ypos(i,j));
    else
        dist = 99999;
        if xpos(i,j)<square_xlo && ypos(i,j)<square_ylo
            dist = -sqrt((xpos(i,j)-square_xlo)^2+(ypos(i,j)-square_ylo)^2);
        elseif xpos(i,j)>=square_xhi && ypos(i,j)<square_ylo
            dist = -sqrt((xpos(i,j)-square_xhi)^2+(ypos(i,j)-square_ylo)^2);
        elseif xpos(i,j)<square_xlo && ypos(i,j)>=square_yhi
            dist = -sqrt((xpos(i,j)-square_xlo)^2+(ypos(i,j)-square_yhi)^2);
        elseif xpos(i,j)>=square_xhi && ypos(i,j)>=square_yhi
            dist = -sqrt((xpos(i,j)-square_xhi)^2+(ypos(i,j)-square_yhi)^2);
        elseif xpos(i,j)<square_xlo
            dist = xpos(i,j)-square_xlo;
        elseif xpos(i,j)>square_xhi
            dist = square_xhi-xpos(i,j);
        elseif ypos(i,j)<square_ylo
            dist = ypos(i,j)-square_ylo;
        elseif ypos(i,j)>square_yhi
            dist = square_yhi-ypos(i,j);
        end
    end
    dext(i,j) = dist;
end
end


if (1)
    figure;
    contour(xpos(lo(1):hi(1),lo(2):hi(2)),ypos(lo(1):hi(1),lo(2):hi(2)), ...
    d0(lo(1):hi(1),lo(2):hi(2)), 'ShowText','on');
    colorbar;
    title('Raw');
    axis equal;
    hold on;
    contour(xpos(lo(1):hi(1),lo(2):hi(2)),ypos(lo(1):hi(1),lo(2):hi(2)), ...
    dext(lo(1):hi(1),lo(2):hi(2)), 'ShowText','on');
    hold off;
end

dd = LSReinit2D(d0, lo,hi,ng,dh);
if (1)
    figure;
    contour(xpos(lo(1):hi(1),lo(2):hi(2)),ypos(lo(1):hi(1),lo(2):hi(2)), ...
    dd(lo(1):hi(1),lo(2):hi(2)), 'ShowText','on');
    colorbar;
    title('Re-initialized');
    axis equal;
    hold on;
    contour(xpos(lo(1):hi(1),lo(2):hi(2)),ypos(lo(1):hi(1),lo(2):hi(2)), ...
    dext(lo(1):hi(1),lo(2):hi(2)), 'ShowText','off');
    hold off;
end


