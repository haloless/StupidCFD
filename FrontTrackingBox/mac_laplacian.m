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

% ## mac_laplacian

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-06-26

function [ lap ] = mac_laplacian (nx,ny, dx,dy, bc)
    
    
    
    % BC is simply ignored
    
    FTGlobals;
    
    N = nx * ny;
    
    % % 5-point discretization
    % lap = spalloc(N,N,N*5);
    
    % Ae = zeros(nx,ny);
    % Ae(1:nx-1,1:ny) = 1/dx^2;
    % Aw = zeros(nx,ny);
    % Aw(2:nx,1:ny) = 1/dx^2;
    % An = zeros(nx,ny);
    % An(1:nx,1:ny-1) = 1/dy^2;
    % As = zeros(nx,ny);
    % As(1:nx,2:ny) = 1/dy^2;   
    
    lap = sparse(N,N); 
    stride = nx;
    
    idx = 0;
    for j = 1:ny
        for i = 1:nx
            idx = idx + 1;
            idw = idx - 1;
            ide = idx + 1;
            ids = idx - stride;
            idn = idx + stride;
            
            rl = EdgeRs(i);
            rr = EdgeRs(i+1);
            rc = 0.5 * (rl + rr);
            
            if (i ~= 1)
                % coef = (rl/rc) / dx^2;
                coef = rl / dx^2;
                lap(idx,idw) = -coef;
                lap(idx,idx) = lap(idx,idx) + coef;
            end
            if (j ~= 1)
                coef = rc / dy^2;
                lap(idx,ids) = -coef;
                lap(idx,idx) = lap(idx,idx) + coef;
            end
            if (j ~= ny)
                coef = rc / dy^2;
                lap(idx,idn) = -coef;
                lap(idx,idx) = lap(idx,idx) + coef;
            end
            if (i ~= nx)
                % coef = (rr/rc) / dx^2;
                coef = rr / dx^2;
                lap(idx,ide) = -coef;
                lap(idx,idx) = lap(idx,idx) + coef;
            end
        end
    end
    
    % set reference point
    lap(1,:) = 0;
    lap(:,1) = 0;
    lap(1,1) = 1;
end


