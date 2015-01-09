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

% ## EBPPERhs

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-29

function [ rhs ] = EBPPERhs (ustar,vstar,nx,ny,dx,dy,dt, bx,by)
% Description

EBGlobals;

I = 2:nx+1;
J = 2:ny+1;
% div_vel = 1/dx * (bx(I+1,J).*ustar(I+1,J) - bx(I,J).*ustar(I,J)) ...
    % + 1/dy * (by(I,J+1).*vstar(I,J+1) - by(I,J).*vstar(I,J));
div_vel = 1/dx*(ustar(I+1,J)-ustar(I,J)) + 1/dy*(vstar(I,J+1)-vstar(I,J));
rhs = -1/dt * rho * div_vel;

% TODO BC correction

% return as a vector
rhs = reshape(rhs,nx*ny,1);
return
end
