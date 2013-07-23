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

% ## WENO_bc

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-16

function [ un ] = WENO_bc (un, nx, ngrow)

% periodic BC
nbc = ngrow + 1;
% high BC
un(nx-nbc+1:nx) = un(nbc+1:nbc+nbc);
% low BC
un(1:nbc) = un(nx-nbc+1-nbc:nx-nbc);

% un(3) = un(nx-2);
% un(nx-2) = un(3);
% un(nx-1:nx) = un(4:5);
% un(1:2) = un(nx-4:nx-3);

% ubc = 0.5*(un(nbc) + un(nx-ngrow));
% un(nbc) = ubc;
% un(nx-ngrow) = ubc;
% un(nx-ngrow+1:nx) = un(nbc+1:nbc+ngrow);
% un(1:ngrow) = un(nx-ngrow-ngrow:nx-ngrow-1);
% un(nx-1:nx) = un(4:5);
% un(1:2) = un(nx-4:nx-3);


return
end



