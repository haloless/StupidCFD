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

function [ nlevel ] = SemiMGInit(nx,ny)

SemiMGGlobals;

% clear SMG
smg_num_level = -1;
smg_levels = [];

%
nlevel = 1;
cdir = [ -1 ];

ncx = nx;
ncy = ny;
while (1)
    ok = 0;
    
    nlenmin = 1;
    
    nc = floor(ncx/2);
    if (nc >= nlenmin)
        nlevel = nlevel + 1;
        cdir(end+1) = 1;
        ncx = floor(ncx/2);
        ok = 1;
    end
    
    nc = floor(ncy/2);
    if (nc >= nlenmin)
        nlevel = nlevel+1;
        cdir(end+1) = 2;
        ncy = floor(ncy/2);
        ok = 1;
    end
    
    if (~ok) 
        break;
    end
end

%
smg_num_level = nlevel;

ncx = nx;
ncy = ny;

%
% for lvl = nlevel:-1:1
for lvl = 1:nlevel
    % idx = nlevel-lvl+1;
    idx = lvl;
    
    if (cdir(idx) == 1)
        ncx = floor(ncx/2);
    elseif (cdir(idx) == 2)
        ncy = floor(ncy/2);
    end
    
    % coarsen direction (-1 on finest)
    smg_levels(lvl).cdir = cdir(idx);
    % cell number 
    smg_levels(lvl).ncell = [ ncx ncy ];
    
    % A matrix
    smg_levels(lvl).acen = zeros(ncx,ncy);
    smg_levels(lvl).axlo = zeros(ncx,ncy);
    smg_levels(lvl).axhi = zeros(ncx,ncy);
    smg_levels(lvl).aylo = zeros(ncx,ncy);
    smg_levels(lvl).ayhi = zeros(ncx,ncy);
    
    % P
    smg_levels(lvl).plo = zeros(ncx,ncy);
    smg_levels(lvl).phi = zeros(ncx,ncy);
    smg_levels(lvl).tcen = zeros(ncx,ncy);
    
    % relaxation
    smg_levels(lvl).relax_weight = 2.0/3.0;
    if (lvl == nlevel)
        smg_levels(lvl).relax_weight = 1.0;
    end
    
    % active
    smg_levels(lvl).active = 1;
    
    % data vector
    smg_levels(lvl).xvec = zeros(ncx,ncy);
    smg_levels(lvl).bvec = zeros(ncx,ncy);
    smg_levels(lvl).rvec = zeros(ncx,ncy);
    smg_levels(lvl).evec = zeros(ncx,ncy);
end

return
end





