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

% ## LSHeaviside

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-04-28

function [ hs ] = LSHeaviside (dist,cutoff)
% Description

if (dist > cutoff)
    hs = 1.0;
elseif (dist < -cutoff)
    hs = 0.0;
else
    phi = dist / cutoff;
    hs = 0.5 * (1.0 + phi + sin(pi*phi)/pi);
end

return
end
