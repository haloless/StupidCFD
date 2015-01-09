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

% ## SumPoly

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-01-18

function [ pw ] = SumPoly (ps, ws)
% Description


if (nargin==2 && size(ps,2)~=length(ws))
    error('SumPoly: weights invalid')
end

if (nargin == 1)
    pw = sum(ps,2);
elseif (nargin == 2)
    if (isrow(ws)); ws = ws'; end
    pw = ps * ws;
end

return
end


