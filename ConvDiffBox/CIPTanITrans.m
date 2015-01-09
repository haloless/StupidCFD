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

% ## CIPTanITrans

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-06

function [ f ] = CIPTanITrans (tf,tan_coef)

if ~exist('tan_coef','var')
    tan_coef = 0.9;
end

f = atan(tf) / (tan_coef*pi) + 0.5;
% f = tanh(tf)/(2*tan_coef) + 0.5;

return
end
