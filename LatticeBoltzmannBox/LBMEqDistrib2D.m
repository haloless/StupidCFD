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

% ## LBMEqDistrib

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-19

function [ fEq ] = LBMEqDistrib2D (rho,u,v,nx,ny, qwgt,qex,qey, qc,cs2)
% Description

cs4 = cs2^2;
nlink = length(qwgt);

fEq = zeros(length(rho),nlink);
vel2 = u.^2 + v.^2;

for i = 1:nlink
    eu = qc * (qex(i)*u + qey(i)*v);
    fEq(:,i) = qwgt(i) * rho .* (1 + 1/cs2*eu + 1/(2*cs4)*(eu.^2) - 1/(2*cs2)*vel2);
end


return
end
