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

% ## RBF_GA

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-28

function [ R,dRdr ] = RBFunc_GA (rs, epsilon)
% Description
% @(r,eps) exp(-(eps.*r).^2);

e2 = epsilon^2;

R = exp(-e2 * rs.^2);
dRdr = -2*e2 * rs .* R;
d2Rdr2 = (-2*e2 + 4*e2^2*rs.^2) .* R;



return
end
