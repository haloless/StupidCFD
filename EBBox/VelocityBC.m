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

% ## VelocityBC

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-07

function [ umac,vmac ] = VelocityBC (umac,vmac,nx,ny)

global UIn;

% umac(1,:) = 0;

% x low
umac(1,2:ny+1) = UIn;
umac(2,2:ny+1) = UIn;
vmac(1,2:ny+2) = -vmac(2,2:ny+2);


% y low
umac(2:nx+2,1) = -umac(2:nx+2,2);
vmac(2:nx+1,1) = vmac(2:nx+1,3);
vmac(2:nx+1,2) = 0;

% y high
umac(2:nx+2,ny+2) = -umac(2:nx+2,ny+1);
vmac(2:nx+1,ny+3) = vmac(2:nx+1,ny+1);
vmac(2:nx+1,ny+2) = 0;

% x high
% Neumann
umac(nx+3,2:ny+1) = umac(nx+2,2:ny+1); 
vmac(nx+2,2:ny+2) = vmac(nx+1,2:ny+2);




return
end
