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

% ## Euler1D_PrimToCons

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-01-03

function [ ucons ] = Euler1D_PrimToCons (uprim, ucons, lo,hi)

Euler1D_globals;

irange = lo:hi;

rho = uprim(QRHO,irange);
vx = uprim(QVX,irange);
p = uprim(QPRES,irange);

ucons(URHO,irange) = rho;
ucons(UMX,irange) = rho .* vx;
ucons(UETOT,irange) = 0.5 .* rho .* vx.^2 + p ./ (GAMMA-1);

return
end
