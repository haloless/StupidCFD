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

% ## RiderKotheVortex

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-10-16

function [ u,v ] = RiderKotheVortex (x,y,t,T)

if (T <= 0)
    coef = 0.0;
else
    coef = cos(pi*t/T);
end

pix = pi * x;
piy = pi * y;

u = -coef .* (sin(pix).^2) .* sin(2*piy);
v =  coef .* (sin(piy).^2) .* sin(2*pix);

return
end
