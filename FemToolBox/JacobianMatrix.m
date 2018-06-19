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

% ## JacobianMatrix

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-30

function [ Jmat,invJ,detJ ] = JacobianMatrix (dN,nodes)
% Description

dNdxi = dN(:,1);
dNdeta = dN(:,2);

xs = nodes(:,1);
ys = nodes(:,2);

Jmat = [dNdxi'*xs, dNdxi'*ys; dNdeta'*xs, dNdeta'*ys];
invJ = inv(Jmat);
detJ = det(Jmat);

return
end
