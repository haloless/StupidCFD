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

% ## ../Utils/easy_vorticity

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-09

function [ vort ] = easy_vorticity (xs,ys,nx,ny,dx,dy,u,v)
% Description
% vorticity = dv/dx - du/dy
% using cell velocity on [-1:nx+1] by [-1:ny+1]


I = 2:nx+1;
J = 2:ny+1;
vort = 1/(2*dx) * (v(I+1,J)-v(I-1,J)) - 1/(2*dy) * (u(I,J+1)-u(I,J-1));

return
end
