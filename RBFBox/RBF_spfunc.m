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

% ## RBF_spfunc

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-24

function [ sp_R,sp_dRdx,sp_dRdy ] = RBF_spfunc (x,y,xs,ys,re,neigh)

[R,dRdx,dRdy] = RBF_func(x,y,xs(neigh),ys(neigh),re);

sz = size(xs);

sp_R = sparse(sz(1),sz(2));
sp_dRdx = sparse(sz(1),sz(2));
sp_dRdy = sparse(sz(1),sz(2));

sp_R(neigh) = R;
sp_dRdx(neigh) = dRdx;
sp_dRdy(neigh) = dRdy;

return
end



