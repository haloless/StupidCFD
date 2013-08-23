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

% ## VelocityDiffusion

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-12

function [ Diffu Diffv ] = VelocityDiffusion (umac,vmac,nx,ny,dx,dy,dt)

EBGlobals;

Diffu = zeros(nx+3,ny+2);
I = 2:nx+2;
J = 2:ny+1;
udiff = 1/dx^2 * (umac(I+1,J) - 2*umac(I,J) + umac(I-1,J)) ...
    + 1/dy^2 * (umac(I,J+1) - 2*umac(I,J) + umac(I,J-1));
Diffu(I,J) = nu * udiff;

Diffv = zeros(nx+2,ny+3);
I = 2:nx+1;
J = 2:ny+2;
vdiff = 1/dx^2 * (vmac(I+1,J) - 2*vmac(I,J) + vmac(I-1,J)) ...
    + 1/dy^2 * (vmac(I,J+1) - 2*vmac(I,J) + vmac(I,J-1));
Diffv(I,J) = nu * vdiff;


return
end
