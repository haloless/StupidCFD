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

function [ ] = SemiMGSetup(Ac,Aw,Ae,As,An)

SemiMGGlobals;

nlevel = smg_num_level;

% set matrix for finest level
smg_levels(1).acen = Ac;
smg_levels(1).axlo = Aw;
smg_levels(1).axhi = Ae;
smg_levels(1).aylo = As;
smg_levels(1).ayhi = An;

% generate P and RAP
for level = 1:nlevel
    if (level > 1)
        SemiMGGenMat(level);
    end
    if (level < nlevel)
        SemiMGGenInterp(level);
    end
end


return
end





