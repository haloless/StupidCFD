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

% ## VOFWLIC_normal2d

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-10-04

function [ normalx,normaly ] = VOFWLIC_normal2d (f,nx,ny,dx,dy)

I = 1:nx+1;
J = 1:ny+1;
nodenx = 1/(2*dx) * (f(I+1,J)-f(I,J) + f(I+1,J+1)-f(I,J+1));
nodeny = 1/(2*dy) * (f(I,J+1)-f(I,J) + f(I+1,J+1)-f(I+1,J));

normalx = zeros(nx+2,ny+2);
normaly = zeros(nx+2,ny+2);

I = 1:nx;
J = 1:ny;
normalx(2:nx+1,2:ny+1) = 0.25 * (nodenx(I,J) + nodenx(I+1,J) + nodenx(I,J+1) + nodenx(I+1,J+1));
normaly(2:nx+1,2:ny+1) = 0.25 * (nodeny(I,J) + nodeny(I+1,J) + nodeny(I,J+1) + nodeny(I+1,J+1));

% normals need a Neumann BC
normalx = VOFWLIC_scalarBC2d(normalx,nx,ny);
normaly = VOFWLIC_scalarBC2d(normaly,nx,ny);


return
end
