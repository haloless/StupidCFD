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

% ## solve_pressure

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-03

function [ unew,vnew,pnew ] = solve_pressure (ustar,vstar,pstar)

simple_globals;

N = nx * ny;
Lp = spalloc(N,N,N*5);
% rhs = zeros(N,1);
rhs = dy * diff(ustar(1:nx+1,2:ny+1)) + dx * diff(vstar(2:nx+1,1:ny+1)')';
rhs = reshape(-rho * rhs, N);

stride = nx;
idx = 0;
for j = 2:ny+1
    for i = 2:nx+1
        idx = idx + 1;
        ide = idx + 1;
        idw = idx - 1;
        idn = idx + stride;
        ids = idx - stride;
        
        % RHS
        % rhs(idx) = -rho*((ustar(i,j)-ustar(i-1,j))*dy+(vstar(i,j)-vstar(i,j-1))*dx);
        
        ape = 0;
        apw = 0;
        apn = 0;
        aps = 0;
        % assemble
        if (j ~= 2)
            aps = rho * dx^2 / AVp(i,j-1);
            Lp(idx,ids) = -aps;
        end
        if (i ~= 2)
            apw = rho * dy^2 / AUp(i-1,j);
            Lp(idx,idw) = -apw;
        end
        if (i ~= nx+1)
            ape = rho * dy^2 / AUp(i,j);
            Lp(idx,ide) = -ape;
        end
        if (j ~= ny+1)
            apn = rho * dx^2 / AVp(i,j);
            Lp(idx,idn) = -apn;
        end
        Lp(idx,idx) = ape + apw + apn + aps;
    end
end

% inject ref. pressure
Lp(1,:) = 0;
Lp(1,1) = 1;
rhs(1) = 0;

sol = Lp \ rhs;

% restore pressure-correction
pdash = zeros(nx+2,ny+2);
pdash(2:nx+1,2:ny+1) = reshape(sol,nx,ny);
pdash = apply_pres_bc(pdash);

% new velocity
udash = zeros(nx+1,ny+2);
udash(2:nx,2:ny+1) = -dy * diff(pdash(2:nx+1,2:ny+1)) ./ AUp(2:nx,2:ny+1);
unew = ustar + udash;
% unew = ustar;
% for j = 2:ny+1
    % for i = 2:nx
        % udash = -dy / AUp(i,j) * (pdash(i+1,j)-pdash(i,j));
        % unew(i,j) = unew(i,j) + udash;
    % end
% end
unew = apply_umac_bc(unew);

vdash = zeros(nx+2,ny+1);
vdash(2:nx+1,2:ny) = -dx * diff(pdash(2:nx+1,2:ny+1)')' ./ AVp(2:nx+1,2:ny);
vnew = vstar + vdash;
% vnew = vstar;
% for j = 2:ny
    % for i = 2:nx+1
        % vdash = -dx / AVp(i,j) * (pdash(i,j+1)-pdash(i,j));
        % vnew(i,j) = vnew(i,j) + vdash;
    % end
% end
vnew = apply_vmac_bc(vnew);

% new pressure
pnew = pstar + prelax * pdash;
pnew = apply_pres_bc(pnew);

return
end
