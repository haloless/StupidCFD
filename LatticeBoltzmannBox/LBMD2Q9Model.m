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

% ## LBMD2Q9Model
% ## 7   3   6
% ##   \ | /
% ## 4 - 1 - 2
% ##   / | \
% ## 8   5   9

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-15

function [ qwgt,qex,qey,qord,qopp ] = LBMD2Q9Model ()

qwgt = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36];
qex  = [0, 1, 0, -1, 0, 1, -1, -1, 1];
qey  = [0, 0, 1, 0, -1, 1, 1, -1, -1];
qord = [1, 2, 3, 4, 5,  6, 7, 8,  9];
qopp = [1, 4, 5, 2, 3,  8, 9, 6,  7];

return
end
