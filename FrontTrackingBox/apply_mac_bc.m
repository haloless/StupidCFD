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

% ## apply_bc

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-06-26

function [ umac, vmac ] = apply_mac_bc (umac,vmac, nx,ny, BC)
    
    FTGlobals;
    
    global UIn
    
    if (1)

    % % x-low
    % umac(1,2:ny+1) = 0;
    % vmac(2:nx+1,ny+1) = 0;
    % % x-high
    % umac(nx+1,2:ny+1) = 0;
    % vmac(nx+2,1:ny+1) = -vmac(nx+1,1:ny+1);
    % % y-low
    % umac(1:nx+1,1) = -umac(1:nx+1,2);
    % vmac(2:nx+1,1) = 0;
    % % y-high
    % umac(1:nx+1,ny+2) = 2*ULid - umac(1:nx+1,ny+1);
    % vmac(1,1:ny+1) = -vmac(2,1:ny+1);
    % elseif CoordSys == 1
        % r-low
        umac(1,2:ny+1) = 0;
        % r-high
        umac(nx+1,2:ny+1) = 0;
        % z-low
        umac(1:nx+1,1) = -umac(1:nx+1,2);
        % z-high
        umac(1:nx+1,ny+2) = umac(1:nx+1,ny+1);
        
        % z-low
        vmac(2:nx+1,1) = UIn;
        % z-high
        vmac(2:nx+1,ny+1) = vmac(2:nx+1,ny);
        if (1)
            % correct outflow
            inflow = sum(vmac(2:nx+1,1) .* CellRs(2:nx+1));
            outflow = sum(vmac(2:nx+1,ny+1) .* CellRs(2:nx+1));
            vmac(2:nx+1,ny+1) = vmac(2:nx+1,ny+1) * (inflow/(outflow+1e-30));
        end
        % r-low
        vmac(1,1:ny+1) = vmac(2,1:ny+1);
        % r-high
        vmac(nx+2,1:ny+1) = -vmac(nx+1,1:ny+1);
        
    end
    
end
