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

% ## VOFWLIC_scalarBC2d

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-10-04

function [ phi ] = VOFWLIC_scalarBC2d (phi,nx,ny)

phi(1,:) = phi(2,:);
phi(nx+2,:) = phi(nx+1,:);

phi(:,1) = phi(:,2);
phi(:,ny+2) = phi(:,ny+1);

return
end
