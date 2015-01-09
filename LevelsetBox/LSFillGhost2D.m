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

% ## LSFillBC

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-04-28

function [ dd ] = LSFillGhost2D (dd,lo,hi,ng)

LSGlobals2D;

if (ng < 1); return; end

% fill homogeneous BC
for ii = 1:ng
    dd(lo(1)-ii,:) = dd(lo(1),:);
    dd(hi(1)+ii,:) = dd(hi(1),:);
end
for jj = 1:ng
    dd(:,lo(2)-jj) = dd(:,lo(2));
    dd(:,hi(2)+jj) = dd(:,hi(2));
end

return
end
