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

function [ fine ] = SemiMGProlong(level, crse)

SemiMGGlobals;

nlevel = smg_num_level;
if (level == nlevel)
    error('SemiMG: Prolong: called on coarsest level');
end

% coarsen direction
cdir = smg_levels(level+1).cdir;
if (cdir<1 || cdir>2)
    error('SemiMG: Prolong: CDIR invalid');
end

% crse cell
ncell = smg_levels(level+1).ncell;
ncx = ncell(1);
ncy = ncell(2);
% fine cell
ncell = smg_levels(level).ncell;
nx = ncell(1);
ny = ncell(2);
% check cell match
if (cdir == 1) 
    if (ny ~= ncy)
        error('SemiMG: Prolong: CDIR=1 NY not match');
    end
elseif (cdir == 2)
    if (nx ~= ncx)
        error('SemiMG: Prolong: CDIR=2 NX not match');
    end
end

%
fplo = smg_levels(level).plo;
fphi = smg_levels(level).phi;

%
fine = zeros(nx,ny);
if (cdir == 1)
    J = 1:ny;
    % 
    for i = 2:2:nx
        ic = i / 2;
        fine(i,J) = crse(ic,J);
    end
    %
    for i = 1:2:nx
        if (i > 1)
            fine(i,J) = fine(i,J) + fplo(i,J).*fine(i-1,J);
        end
        if (i < nx)
            fine(i,J) = fine(i,J) + fphi(i,J).*fine(i+1,J);
        end
    end
elseif (cdir == 2)
    I = 1:nx;
    %
    for j = 2:2:ny
        jc = j / 2;
        fine(I,j) = crse(I,jc);
    end
    %
    for j = 1:2:ny
        if (j > 1)
            fine(I,j) = fine(I,j) + fplo(I,j).*fine(I,j-1);
        end
        if (j < ny)
            fine(I,j) = fine(I,j) + fphi(I,j).*fine(I,j+1);
        end
    end
end


return
end





