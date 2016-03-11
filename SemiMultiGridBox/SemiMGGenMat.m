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

function [ ] = SemiMGGenMat(level)

SemiMGGlobals;

nlevel = smg_num_level;
if (level == 1)
    error('SemiMG: RAP: called on finest level');
end

% coarsen direction
cdir = smg_levels(level).cdir;
if (cdir<1 || cdir>2)
    error('SemiMG: RAP: CDIR invalid');
end

% crse cell
ncell = smg_levels(level).ncell;
ncx = ncell(1);
ncy = ncell(2);
% fine cell
ncell = smg_levels(level-1).ncell;
nx = ncell(1);
ny = ncell(2);
% check
if (cdir == 1)
    if (ncy ~= ny)
        error('SemiMG: RAP: CDIR=1 NY not match');
    end
elseif (cdir == 2)
    if (ncx ~= nx)
        error('SemiMG: RAP: CDIR=2 NX not match');
    end
end


% fine entries
facen = smg_levels(level-1).acen;
faxlo = smg_levels(level-1).axlo;
faxhi = smg_levels(level-1).axhi;
faylo = smg_levels(level-1).aylo;
fayhi = smg_levels(level-1).ayhi;
%
ftcen = smg_levels(level-1).tcen;
fplo = smg_levels(level-1).plo;
fphi = smg_levels(level-1).phi;

% crse entries
acen = zeros(ncx,ncy);
axlo = zeros(ncx,ncy);
axhi = zeros(ncx,ncy);
aylo = zeros(ncx,ncy);
ayhi = zeros(ncx,ncy);

if (cdir == 1)
    for jc = 1:ncy
    for ic = 1:ncx
        i = ic * 2;
        j = jc;
        
        % coef. in cdir
        if (i > 1)
            aw = faxlo(i,j) * fplo(i-1,j);
        else
            aw = 0;
        end
        if (i < nx)
            ae = faxhi(i,j) * fphi(i+1,j);
        else
            ae = 0;
        end
        
        % coef. other
        as = faylo(i,j);
        an = fayhi(i,j);
        if (i > 1)
            as = as + 0.5*faylo(i-1,j);
            an = an + 0.5*fayhi(i-1,j);
        end
        if (i < nx)
            as = as + 0.5*faylo(i+1,j);
            an = an + 0.5*fayhi(i+1,j);
        end
        
        % center
        ac = ftcen(i,j);
        if (i > 1)
            ac = ac + faxlo(i,j) * fphi(i-1,j);
        end
        if (i < nx)
            ac = ac + faxhi(i,j) * fplo(i+1,j);
        end
        ac = ac - (as+an);
        
        acen(ic,jc) = ac;
        axlo(ic,jc) = aw;
        axhi(ic,jc) = ae;
        aylo(ic,jc) = as;
        ayhi(ic,jc) = an;
    end
    end
elseif (cdir == 2)
    for jc = 1:ncy
    for ic = 1:ncx
        i = ic;
        j = jc * 2;
        
        % coef. in cdir
        if (j > 1)
            as = faylo(i,j) * fplo(i,j-1);
        else
            as = 0;
        end
        if (j < ny)
            an = fayhi(i,j) * fphi(i,j+1);
        else
            an = 0;
        end
        
        % other
        aw = faxlo(i,j);
        ae = faxhi(i,j);
        if (j > 1)
            aw = aw + 0.5*faxlo(i,j-1);
            ae = ae + 0.5*faxhi(i,j-1);
        end
        if (j < ny)
            aw = aw + 0.5*faxlo(i,j+1);
            ae = ae + 0.5*faxhi(i,j+1);
        end
        
        %
        ac = ftcen(i,j);
        if (j > 1)
            ac = ac + faylo(i,j) * fphi(i,j-1);
        end
        if (j < ny)
            ac = ac + fayhi(i,j) * fplo(i,j+1);
        end
        ac = ac - (aw+ae);
        
        acen(ic,jc) = ac;
        axlo(ic,jc) = aw;
        axhi(ic,jc) = ae;
        aylo(ic,jc) = as;
        ayhi(ic,jc) = an;
    end
    end
end

% 
smg_levels(level).acen = acen;
smg_levels(level).axlo = axlo;
smg_levels(level).axhi = axhi;
smg_levels(level).aylo = aylo;
smg_levels(level).ayhi = ayhi;


return
end





