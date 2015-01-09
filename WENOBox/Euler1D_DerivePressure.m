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

% ## Euler1D_DerivePressure

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-01-06

function [ pres ] = Euler1D_DerivePressure (ucons)

Euler1D_globals;

ek = 0.5 * ucons(UMX)^2 / ucons(URHO);
pres = (GAMMA-1) * (ucons(UETOT) - ek);

return
end
