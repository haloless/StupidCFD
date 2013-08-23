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

% ## MonomialBasis1D

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-29

function [ p,dpdx ] = MonomialBasis1D (p_order,xs)

if iscolumn(xs); xs=xs'; end;

n = length(xs);
e0 = zeros(1,n);
e1 = ones(1,n);

switch p_order
case {0}
    p = e1;
    dpdx = e0;
case {1}
    p = [e1; xs];
    dpdx = [e0; e1];
case {2}
    p = [e1; xs; xs.^2];
    dpdx = [e0; e1; 2*xs];
case {3}
    p = [e1; xs; xs.^2; xs.^3];
    dpdx = [e0; e1; 2*xs; 3*xs.^2];
otherwise
    error('Unsupported basis order=%d',p_order);
end


return
end




