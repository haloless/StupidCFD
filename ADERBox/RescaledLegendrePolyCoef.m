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

% ## RescaledLegendrePolyCoef

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-01-20

function [ coef ] = RescaledLegendrePolyCoef (n)

switch n
case 0
    coef = [1];
case 1
    % psi = 2*xi - 1;
    coef = [2, -1];
case 2
    % psi = 6*xi.^2 - 6*xi + 1;
    coef = [6, -6, 1];
case 3
    % psi = 20*xi.^3 -30*xi.^2 + 12*xi - 1;
    coef = [20, -30, 12, -1];
case 4
    % psi = 70*xi.^4 - 140*xi.^3 + 90*xi.^2 - 20*xi + 1;
    coef = [70, -140, 90, -20, 1];
case 5
    % psi = 252*xi.^5 - 630*xi.^4 + 560*xi.^3 - 210*xi.^2 + 30*xi - 1;
    coef = [252, -630, 560, -210, 30, -1];
otherwise
    error('RescaledLegendrePolyVal: invalid N');
end

coef = coef';

return
end
