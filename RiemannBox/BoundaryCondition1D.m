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

% ## BoundaryCondition1D

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-01

function [ u ] = BoundaryCondition1D (u,ncell,ngrow,bc_type)
% Description

switch bc_type
case 0 % periodic BC
    u(1:ngrow) = u(ncell+1:ncell+ngrow);
    u(ncell+ngrow+1:ncell+ngrow*2) = u(ngrow+1:ngrow*2);
otherwise
    error('Unsupported BC type: %d',bc_type);
end

return
end
