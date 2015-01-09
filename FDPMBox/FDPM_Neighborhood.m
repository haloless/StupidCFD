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

% ## FDPM_Neighborhood

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-20

function [ neigh,dist,rx,ry,re ] = FDPM_Neighborhood (target,nodes,re)
% Description
% If $TARGET is an integer, then used as an index in $NODES
% otherwise used as coordinate

xc = nodes(target,1);
yc = nodes(target,2);

rx = nodes(:,1)-xc;
ry = nodes(:,2)-yc;
dist = sqrt(rx.^2 + ry.^2);

if (~isscalar(re))
    % re is a vector, to be symmetrized
    rc = re(target);
    re = 1/2 * (re + rc);
end

flag = dist < re;
% exclude target particle itself
flag(target) = 0;
neigh = find(flag);

return
end
