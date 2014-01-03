% ## Copyright (C) 2014 homu
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

% ## Euler1D_ConsToPrim

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-01-03

function [ uprim ] = Euler1D_ConsToPrim (ucons, uprim, lo,hi)

% cons = [rho, mx, E]
% prim = [rho, vx, p]

Euler1D_globals;

irange = lo:hi;

rho = ucons(URHO,irange);
vx = ucons(UMX,irange) ./ rho;
ek = 0.5 .* rho .* vx.^2;

uprim(QRHO,irange) = rho;
uprim(QVX,irange) = vx;
uprim(QPRES,irange) = (GAMMA-1) .* (ucons(3,irange) - ek);

return
end



