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

function [ sol ] = SemiMGRelax(level, sol, rhs)

SemiMGGlobals;

nlevel = smg_num_level;

if (level < nlevel)
    % pre-smooth
    for ismooth = 1:smg_nu1
        sol = SemiMGSmooth(level, sol, rhs);
    end
    
    resid = SemiMGResid(level, sol, rhs);
    
    % restriction
    bcrse = SemiMGRestrict(level+1, resid);
    xcrse = zeros(smg_levels(level+1).ncell);
    
    for irelax = 1:smg_nu0
        xcrse = SemiMGRelax(level+1, xcrse, bcrse);
    end
    
    % prolongation
    corr = SemiMGProlong(level, xcrse);
    sol = sol + corr;
    
    % post-smooth
    for ismooth = 1:smg_nu2
        sol = SemiMGSmooth(level, sol, rhs);
    end
else
    % bottom
    sol = SemiMGSmooth(level, sol, rhs);
end



return
end





