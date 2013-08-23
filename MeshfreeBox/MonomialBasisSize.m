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

% ## MonomialBasisSize

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-29

function [ m ] = MonomialBasisSize (ndim,porder)
% Description

switch ndim
case {1}
    m = porder + 1;
case {2}
    m = (porder+1)*(porder+2) / 2;
case {3}
    m = (porder+1)*(porder+2)*(porder+3) / 6;
otherwise
    error('Ndim must <= 3');
end

return
end
