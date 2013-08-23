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

% ## VelocityPredictor

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-07

function [ Hx Hy ] = VelocityPredictor (umac,vmac,nx,ny,dx,dy,dt)

EBGlobals;

% x direction
Hx = zeros(nx+3,ny+2);
I = 2:nx+2;
J = 2:ny+1;

ue = 0.5 * (umac(I+1,J) + umac(I,J));
uw = 0.5 * (umac(I-1,J) + umac(I,J));
un = 0.5 * (umac(I,J+1) + umac(I,J));
us = 0.5 * (umac(I,J-1) + umac(I,J));
vn = 0.5 * (vmac(I-1,J+1) + vmac(I,J+1));
vs = 0.5 * (vmac(I-1,J) + vmac(I,J));

uadv = 1/dx * (max(ue,0).*umac(I,J) + min(ue,0).*umac(I+1,J)) ...
    - 1/dx * (max(uw,0).*umac(I-1,J) + min(uw,0).*umac(I,J)) ...
    + 1/dy * (max(vn,0).*umac(I,J) + min(vn,0).*umac(I,J+1)) ...
    - 1/dy * (max(vs,0).*umac(I,J-1) + min(vs,0).*umac(I,J));

udiff = 1/dx^2 * (umac(I+1,J) - 2*umac(I,J) + umac(I-1,J)) ...
    + 1/dy^2 * (umac(I,J+1) - 2*umac(I,J) + umac(I,J-1));

Hx(I,J) = -uadv + nu*udiff;

% clear ue uw un us vn vs;

% y direction
Hy = zeros(nx+2,ny+3);
I = 2:nx+1;
J = 2:ny+2;

ue = 0.5 * (umac(I+1,J-1) + umac(I+1,J));
uw = 0.5 * (umac(I,J-1) + umac(I,J));
ve = 0.5 * (vmac(I+1,J) + vmac(I,J));
vw = 0.5 * (vmac(I-1,J) + vmac(I,J));
vn = 0.5 * (vmac(I,J+1) + vmac(I,J));
vs = 0.5 * (vmac(I,J-1) + vmac(I,J));

vadv = 1/dx * (max(ue,0).*vmac(I,J) + min(ue,0).*vmac(I+1,J)) ...
    - 1/dx * (max(uw,0).*vmac(I-1,J) + min(uw,0).*vmac(I,J)) ...
    + 1/dy * (max(vn,0).*vmac(I,J) + min(vn,0).*vmac(I,J+1)) ...
    - 1/dy * (max(vs,0).*vmac(I,J-1) + min(vs,0).*vmac(I,J));

vdiff = 1/dx^2 * (vmac(I+1,J) - 2*vmac(I,J) + vmac(I-1,J)) ...
    + 1/dy^2 * (vmac(I,J+1) - 2*vmac(I,J) + vmac(I,J-1));

Hy(I,J) = -vadv + nu*vdiff;


return
end

