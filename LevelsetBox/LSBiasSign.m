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

% ## LSBiasSign

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-04-28

function [ bias ] = LSBiasSign (dist,dh)

LSGlobals2D;

if (1)
hs = LSHeaviside(dist, dh*mass_spread);
else
hs = dist/(dh*0.5);
end
bias = 2.0 * (hs - 0.5);

return
end
