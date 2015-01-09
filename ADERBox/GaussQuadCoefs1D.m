% ## Copyright (C) 2014 homu
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

% ## GaussQuadCoefs1D

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-01-19

function [ eta, wgt ] = GaussQuadCoefs1D (np, xlo, xhi)

% points and weights on [-1,+1]
switch np
case {3}
    eta = [-sqrt(3/5); 0; sqrt(3/5)];
    wgt = [5/9; 8/9; 5/9];
case {4}
    eta = [-sqrt((3+2*sqrt(6/5))/7), -sqrt((3-2*sqrt(6/5))/7), ...
        sqrt((3-2*sqrt(6/5))/7), sqrt((3+2*sqrt(6/5))/7)]';
    wgt = [(18-sqrt(30))/36, (18+sqrt(30))/36, ...
        (18+sqrt(30))/36, (18-sqrt(30))/36]';
case {5}
    eta = [-sqrt(5+2*sqrt(10/7))/3, -sqrt(5-2*sqrt(10/7))/3, ...
        0, sqrt(5-2*sqrt(10/7))/3, sqrt(5+2*sqrt(10/7))/3]';
    wgt = [(322-13*sqrt(70))/900, (322+13*sqrt(70))/900, ...
        128/225, (322+13*sqrt(70))/900, (322-13*sqrt(70))/900]';
otherwise
    error('GaussQuadCoefs1D: Np invalid')
end

% shift to expected range
if (nargin == 3)
    eta = xlo + 0.5*(xhi-xlo) .* (eta+1);
    wgt = 0.5*(xhi-xlo) .* wgt;
end

return
end
