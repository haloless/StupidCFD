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

% ## mac_rhs

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-06-26

function [ rhs ] = mac_rhs (ustar,vstar, nx,ny, dx,dy,dt)

FTGlobals;

N = nx * ny;
% rhs = zeros(N,1);

% idx = 0;
% for j = 2:ny+1
    % for i = 2:nx+1
        % idx = idx + 1;
        
        % divu = (ustar(i,j)-ustar(i-1,j)) / dx + (vstar(i,j)-vstar(i,j-1)) / dy;
        % rhs(idx) = -divu / dt;
    % end
% end


% rhs = 1/dx * (ustar(2:nx+1,2:ny+1) - ustar(1:nx,2:ny+1)) + ... 
    % 1/dy* (vstar(2:nx+1,2:ny+1) - vstar(2:nx+1,1:ny));

rf = repmat(EdgeRs, [1,ny]);
rc = repmat(CellRs(2:nx+1), [1,ny]);

rhs = 1/dx * (rf(2:nx+1,:).*ustar(2:nx+1,2:ny+1) - rf(1:nx,:).*ustar(1:nx,2:ny+1)) ./ rc + ...
    1/dy * (vstar(2:nx+1,2:ny+1) - vstar(2:nx+1,1:ny));
rhs = rhs .* rc;

rhs = -1/dt * rhs;
rhs = reshape(rhs, N, []);

% inject reference pressure
rhs(1) = 0;


return
end



