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

% ## apply_pres_bc

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-06-27

function [ p ] = apply_pres_bc (p, nx,ny, bc)
    % apply pressure BC
    % p(1,2:ny+1) = p(2,2:ny+1);
    % p(nx+2,2:ny+1) = p(nx+1,2:ny+1);
    % p(2:nx+1,1) = p(2:nx+1,2);
    % p(2:nx+1,ny+2) = p(2:nx+1,ny+1);
    
	p(1,:) = p(nx+1,:);
	p(nx+2,:) = p(2,:);
	p(:,1) = p(:,ny+1);
	p(:,ny+2) = p(:,2);
	
    return
end
