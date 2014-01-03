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

% ## VelocityLapRhs

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-26

function [ uRhs,vRhs ] = VelocityLapRhs (ustar,vstar,ubc,vbc,nx,ny,dx,dy,dt)

%
EBGlobals;

% uRhs = zeros(nx+1,ny);
% vRhs = zeros(nx,ny+1);

coefx = nu*dt/2 * 1/dx^2;
coefy = nu*dt/2 * 1/dy^2;

uRhs = ustar(2:nx+2,2:ny+1);
% x-low
uxlo = ubc(2,2:ny+1);
uRhs(1,:) = uxlo;
uRhs(2,:) = uRhs(2,:) + coefx*uxlo;
%
uRhs = reshape(uRhs, (nx+1)*ny,1);

vRhs = vstar(2:nx+1,2:ny+2);
% y-low
vylo = vbc(2:nx+1,2);
vRhs(:,1) = vylo;
vRhs(:,2) = vRhs(:,2) + coefy*vylo;
% y-high
vyhi = vbc(2:nx+1,ny+2);
vRhs(:,ny+1) = vyhi;
vRhs(:,ny) = vRhs(:,ny) + coefy*vyhi;
%
vRhs = reshape(vRhs, nx*(ny+1),1);


return
end
