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

% ## apply_vmac_bc

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-02

function [ v ] = apply_vmac_bc (v)

simple_globals;

% y-low
v(2:nx+1,1) = 0;
% y-high
v(2:nx+1,ny+1) = 0;

% x-low
v(1,1:ny+1) = -v(2,1:ny+1);
% x-high
v(nx+2,1:ny+1) = -v(nx+1,1:ny+1);

return
end
