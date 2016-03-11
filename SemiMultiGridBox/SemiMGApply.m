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

% ## MGGlobals

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-01-25

function [ out ] = SemiMGApply(level, in)

SemiMGGlobals;

ncell = smg_levels(level).ncell;
nx = ncell(1);
ny = ncell(2);

Ac = smg_levels(level).acen;
Aw = smg_levels(level).axlo;
Ae = smg_levels(level).axhi;
As = smg_levels(level).aylo;
An = smg_levels(level).ayhi;

% out = zeros(nx,ny);
% for j = 1:ny
% for i = 1:nx
    % val = Ac(i,j) * in(i,j);
    % if (i > 1)
        % val = val + Aw(i,j)*in(i-1,j);
    % end
    % if (i < nx)
        % val = val + Ae(i,j)*in(i+1,j);
    % end
    % if (j > 1)
        % val = val + As(i,j)*in(i,j-1);
    % end
    % out(i,j) =  + 
% end
% end

out = Ac .* in;

I = (2:nx);
J = (1:ny);
tmp = Aw(I,J) .* in(I-1,J);
out(I,J) = out(I,J) + tmp;

I = (1:nx-1);
J = (1:ny);
tmp = Ae(I,J) .* in(I+1,J);
out(I,J) = out(I,J) + tmp;

I = (1:nx);
J = (2:ny);
tmp = As(I,J) .* in(I,J-1);
out(I,J) = out(I,J) + tmp;

I = (1:nx);
J = (1:ny-1);
tmp = An(I,J) .* in(I,J+1);
out(I,J) = out(I,J) + tmp;

return
end





