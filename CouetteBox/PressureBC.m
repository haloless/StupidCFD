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

% ## PressureBC

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-07

function [ p ] = PressureBC (p,nx,ny)

EBGlobals;

% x low
p(1,2:ny+1) = p(2,2:ny+1);

% x high
% Dirichlet
% p(nx+2,2:ny+1) = POut;
p(nx+2,2:ny+1) = p(nx+1,2:ny+1);

% y low
p(2:nx+1,1) = p(2:nx+1,2);

% y high
p(2:nx+1,ny+2) = p(2:nx+1,ny+1);


return
end
