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

% ## easy_streamline

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-06-30

function [ psi ] = easy_streamfunc (xs, ys, nx, ny, dx, dy, u, v)

% Description:
% u = d(psi)/dy
% v = -d(psi)/dx

dpsidx = -v;
dpsidy = u;

psi = zeros(nx, ny);
for i = 2:nx
    psi(i,1:ny) = psi(i-1,1:ny) - dx * v(i,1:ny);
end
for j = 2:ny
    psi(1:nx,j) = psi(1:nx,j-1) + dy * u(1:nx,j);
end

psi = psi * 0.5;

end
