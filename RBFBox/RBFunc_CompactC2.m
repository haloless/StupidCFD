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

% ## RBFunc_CompactC2

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-29

function [ R,dRdr ] = RBFunc_CompactC2 (r,supp)
% Description

q = r ./ supp;
q2 = q.^2;
q3 = q.^3;
q4 = q.^4;

Rq = (1-q).^5 .* (8 + 40*q + 48*q2 + 25*q3 + 5*q4);
dRdq = -5*(1-q).^4 .* (8 + 40*q + 48*q2 + 25*q3 + 5*q4) + ...
    (1-q).^5 .* (40 + 96*q + 75*q2 + 20*q3);

R = (q<=1) .* Rq;
dRdr = (q<=1) .* dRdq ./ supp;

end
