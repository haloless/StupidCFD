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

% ## WendlandCSRBF

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-31

function [ R dRdr d2Rdr2 ] = WendlandCSRBF (r,supp)
% Description

% C0
phi_2_0 = @(q) (1-q).^2;
% C2
phi_3_1 = @(q) (1-q).^4 .* (4*q+1);
% C4
phi_4_2 = @(q) (1-q).^6 .* (35*q.^2 + 18*q + 3);
% C6
phi_5_3 = @(q) (1-q).^8 .* (32*q.^3 + 25*q.^2 + 8*r + 1);

% smooth = 6;

% switch smooth
% case
% end

end




