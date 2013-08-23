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

% ## BurgersExactRiemann1D

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-01

function [ ustar ] = BurgersExactRiemann1D (ul,ur)
% Description

if (ul > ur)
    % shock wave
    % compute shock speed
    s = 0.5 * (ul + ur);
    if s >= 0
        ustar = ul;
    else
        ustar = ur;
    end
else
    % rarefaction wave
    
    if ul >= 0
        % right supersonic rarefaction
        ustar = ul;
    elseif ur <= 0
        % left supersonic rarefaction
        ustar = ur;
    elseif ul<0 && ur>0
        % transsonic rarefaction
        ustar = 0;
    end
end


return
end
