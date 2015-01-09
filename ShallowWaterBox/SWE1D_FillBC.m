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

% ## SWE1D_FillBC

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-02-04

function [ ucons ] = SWE1D_FillBC (ucons, lo,hi,ng)

SWE1D_Globals;

for comp = 1:NUCONS
    switch bctype(1,comp)
    case {BC_NEU}
        ucons(comp,lo-ng:lo-1) = ucons(comp,lo);
    case {BC_DIR}
        ucons(comp,lo-ng:lo-1) = bcfill(1,comp);
    otherwise
        error('BC type invalid');
    end
    switch bctype(2,comp)
    case {BC_NEU}
        ucons(comp,hi+1:hi+ng) = ucons(comp,hi);
    case {BC_DIR}
        ucons(comp,hi+1:hi+ng) = bcfill(2,comp);
    otherwise
        error('BC type invalid');
    end
end



return
end
