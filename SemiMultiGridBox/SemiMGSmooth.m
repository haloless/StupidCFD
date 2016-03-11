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

% ## MGGlobals

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-01-25

function [ out ] = SemiMGSmooth(level, in, rhs)

SemiMGGlobals;

Ac = smg_levels(level).acen;

% 
w = smg_levels(level).relax_weight;

out = SemiMGApply(level, in);
out = rhs - out;
out = out ./ Ac;

out = in + w.*out;

return
end





