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

% ## apply_umac_bc

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-02

function [ u ] = apply_umac_bc (u)

simple_globals;

% x-low
u(1,2:ny+1) = 0;
% x-high
u(nx+1,2:ny+1) = 0;

% y-low
u(1:nx+1,1) = -u(1:nx+1,2);
% y-high
u(1:nx+1,ny+2) = 2*ULid - u(1:nx+1,ny+1);

return
end




