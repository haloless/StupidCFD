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

% ## LBMStream2d

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-15

function [ fs ] = LBMStream2d (f,nx,ny,qex,qey)

[mx,my] = size(f);

fs = reshape(f,nx,ny);

if (0)
    if (qex == 1)
        tmp = fs(nx,:);
        fs(2:nx,:) = fs(1:nx-1,:);
        fs(1,:) = tmp;
    elseif (qex == -1)
        tmp = fs(1,:);
        fs(1:nx-1,:) = fs(2:nx,:);
        fs(nx,:) = tmp;
    end
    if (qey == 1)
        tmp = fs(:,ny);
        fs(:,2:ny) = fs(:,1:ny-1);
        fs(:,1) = tmp;
    elseif (qey == -1)
        tmp = fs(:,1);
        fs(:,1:ny-1) = fs(:,2:ny);
        fs(:,ny) = tmp;
    end
else
    fs = circshift(fs, [qex qey]);
end

if (mx~=nx || my~=ny)
    fs = reshape(fs, mx,my);
end

return
end
