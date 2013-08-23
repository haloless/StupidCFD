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

% ## GaussPoint

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-30

function [ point,weight ] = GaussPoint (numNode,numGauss)

switch numGauss
    case {4}
        a = 1 / sqrt(3);
        xi = [-a; a; -a; a];
        eta = [-a; -a; a; a];
        point = [xi, eta];
        weight = ones(4,1);
    otherwise
        error('Unsupported Node=%d Gauss=%d',numNode,numGauss);
end



return
end
