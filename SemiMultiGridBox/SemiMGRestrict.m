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

function [ crse ] = SemiMGRestrict(level, fine)

SemiMGGlobals;

nlevel = smg_num_level;
if (level == 1)
    error('SemiMG: Restrict: called on finest level');
end

% crse cell
ncell = smg_levels(level).ncell;
ncx = ncell(1);
ncy = ncell(2);

% coarsen direction
cdir = smg_levels(level).cdir;

% fine cell
ncell = smg_levels(level-1).ncell;
nx = ncell(1);
ny = ncell(2);

%
fplo = smg_levels(level-1).plo;
fphi = smg_levels(level-1).phi;

crse = zeros(ncx,ncy);
if (cdir == 1)
    if (ncy ~= ny)
        error('cdir=1 NY not match');
    end
    
    J = 1:ny;
    for ic = 1:ncx
        i = ic * 2;
        tmp = fine(i,J);
        if (i-1 >= 1)
            tmp = tmp + fphi(i-1,J).*fine(i-1,J);
        end
        if (i+1 <= nx)
            tmp = tmp + fplo(i+1,J).*fine(i+1,J);
        end
        crse(ic,J) = tmp;
    end
elseif (cdir == 2)
    if (ncx ~= nx)
        error('cdir=2 NX not match');
    end
    
    I = 1:nx;
    for jc = 1:ncy
        j = jc * 2;
        tmp = fine(I,j);
        if (j-1 >= 1)
            tmp = tmp + fphi(I,j-1).*fine(I,j-1);
        end
        if (j+1 <= ny)
            tmp = tmp + fplo(I,j+1).*fine(I,j+1);
        end
        crse(I,jc) = tmp;
    end
end


return
end





