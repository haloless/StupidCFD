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

% ## SWE1D_PrimToCons

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-02-02

function [ ucons ] = SWE1D_PrimToCons (uprim)

SWE1D_Globals;

ucons = zeros(NUCONS,size(uprim,2));
%
hs = uprim(QH,:);
us = uprim(QVX,:);
bs = uprim(QBOT,:);
%
ucons(UH,:) = hs;
ucons(UHU,:) = hs .* us;
ucons(UBOT,:) = bs;

return
end
