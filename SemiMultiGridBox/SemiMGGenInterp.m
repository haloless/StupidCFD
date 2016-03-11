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

function [ ] = SemiMGGenInterp(level)

SemiMGGlobals;

nlevel = smg_num_level;
if (level == nlevel)
    error('SemiMG: P: called on coarsest level');
end

%
ncell = smg_levels(level).ncell;
nx = ncell(1);
ny = ncell(2);

% coarsen direction
cdir = smg_levels(level+1).cdir;
if (cdir<1 || cdir>2)
    error('SemiMG: P: CDIR invalid');
end

% matrix entries
acen = smg_levels(level).acen;
axlo = smg_levels(level).axlo;
axhi = smg_levels(level).axhi;
aylo = smg_levels(level).aylo;
ayhi = smg_levels(level).ayhi;

%
tcen = zeros(nx,ny);
if (cdir == 1)
    tcen = acen + aylo + ayhi;
elseif (cdir == 2)
    tcen = acen + axlo + axhi;
end

%
plo = zeros(nx,ny);
phi = zeros(nx,ny);
if (cdir == 1)
    plo = -axlo ./ tcen;
    phi = -axhi ./ tcen;
elseif (cdir == 2)
    plo = -aylo ./ tcen;
    phi = -ayhi ./ tcen;
end

%
smg_levels(level).plo = plo;
smg_levels(level).phi = phi;
smg_levels(level).tcen = tcen;


return
end





