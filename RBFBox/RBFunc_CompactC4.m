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

% ## RBFunc_CompactC4

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-29

function [ R,dRdr ] = RBFunc_CompactC4 (r,supp)

q = r ./ supp;
q2 = q.^2;
q3 = q.^3;
q4 = q.^4;
q5 = q.^5;

Rq = (1-q).^6 .* (6 + 36*q + 82*q2 + 72*q3 + 30*q4 + 5*q5);
dRdq = -6*(1-q).^5 .* (6 + 36*q + 82*q2 + 72*q3 + 30*q4 + 5*q5) + ...
    (1-q).^6 .* (36 + 164*q + 216*q2 + 120*q3 + 25*q4);

R = (q<=1) .* Rq;
dRdr = (q<=1) .* dRdq ./ supp;

return
end


