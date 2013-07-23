% ## Copyright (C) 2013 admin_2
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

% ## MLS_neigh

% ## Author: admin_2 <admin_2@KOSHIZUKA>
% ## Created: 2013-07-19

function [ neigh ] = MLS_neigh (x,y,xs,ys,re)

rs2 = (x-xs).^2 + (y-ys).^2;
if iscolumn(rs2)
	rs2 = rs2';
end

neigh = sparse(rs2<=re^2);

return
end
