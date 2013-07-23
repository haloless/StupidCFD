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

% ## solve_vstar

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-03

function [ vstar ] = solve_vstar (umac,vmac,uold,vold,pstar)

simple_globals;

relax = vrelax;

N = nx * (ny-1);
% H-operator
H = spalloc(N,N,5*N);
% rhs = zeros(N,1);

rhs = zeros(nx+2,ny+1);

% vectorized version
% diffusion coef.
dve = visc * dy / dx;
dvw = visc * dy / dx;
dvn = visc * dx / dy;
dvs = visc * dx / dy;

% advection
fve = rho * dy * 0.5*(uold(2:nx+1,2:ny)+uold(2:nx+1,3:ny+1));
fvw = rho * dy * 0.5*(uold(1:nx,2:ny)+uold(1:nx,3:ny+1));
fvn = rho * dx * 0.5*(vold(2:nx+1,2:ny)+vold(2:nx+1,3:ny+1));
fvs = rho * dx * 0.5*(vold(2:nx+1,2:ny)+vold(2:nx+1,1:ny-1));

ave = zeros(nx+2,ny+1);
avw = zeros(nx+2,ny+1);
avn = zeros(nx+2,ny+1);
avs = zeros(nx+2,ny+1);
ave(2:nx+1,2:ny) = dve + max(-fve,0);
avw(2:nx+1,2:ny) = dvw + max(fvw,0);
avn(2:nx+1,2:ny) = dvn + max(-fvn,0);
avs(2:nx+1,2:ny) = dvs + max(fvs,0);
avp = ave+avw+avn+avs + rho*dx*dy/dt;
% save AVp
AVp(2:nx+1,2:ny) = avp(2:nx+1,2:ny);

% RHS
rhs(2:nx+1,2:ny) = -diff(pstar(2:nx+1,2:ny+1)')' * dx + ...
    rho*dx*dy/dt * vold(2:nx+1,2:ny);
% E
rhs(nx+1,2:ny) = rhs(nx+1,2:ny) + ave(nx+1,2:ny) .* vmac(nx+2,2:ny);
ave(nx+1,2:ny) = 0;
% W
rhs(2,2:ny) = rhs(2,2:ny) + avw(2,2:ny) .* vmac(1,2:ny);
avw(2,2:ny) = 0;
% N
rhs(2:nx+1,ny) = rhs(2:nx+1,ny) + avn(2:nx+1,ny) .* vmac(2:nx+1,ny+1);
avn(2:nx+1,ny) = 0;
% S
rhs(2:nx+1,2) = rhs(2:nx+1,2) + avs(2:nx+1,2) .* vmac(2:nx+1,1);
avs(2:nx+1,2) = 0;
%
rhs(2:nx+1,2:ny) = rhs(2:nx+1,2:ny) + ...
    (1-relax)/relax * avp(2:nx+1,2:ny) .* vold(2:nx+1,2:ny);
rhs = reshape(rhs(2:nx+1,2:ny), N);

idx = 0;
stride = nx;
for j = 2:ny
    for i = 2:nx+1
        idx = idx + 1;
        ide = idx + 1;
        idw = idx - 1;
        idn = idx + stride;
        ids = idx - stride;
        
        if (j ~= 2)
            H(idx,ids) = -avs(i,j);
        end
        if (i ~= 2)
            H(idx,idw) = -avw(i,j);
        end
        H(idx,idx) = avp(i,j) / relax;
        if (i ~= nx+1)
            H(idx,ide) = -ave(i,j);
        end
        if (j ~= ny)
            H(idx,idn) = -avn(i,j);
        end
    end
end


% idx = 0;
% stride = nx;
% for j = 2:ny
    % for i = 2:nx+1
        % idx = idx + 1;
        % ide = idx + 1;
        % idw = idx - 1;
        % idn = idx + stride;
        % ids = idx - stride;
        
        % % diffusion coef.
        % dve = visc * dy / dx;
        % dvw = visc * dy / dx;
        % dvn = visc * dx / dy;
        % dvs = visc * dx / dy;
        
        % % advection coef.
        % fve = rho * 0.5*(uold(i,j)+uold(i,j+1)) * dy;
        % fvw = rho * 0.5*(uold(i-1,j)+uold(i-1,j+1)) * dy;
        % fvn = rho * 0.5*(vold(i,j)+vold(i,j+1)) * dx;
        % fvs = rho * 0.5*(vold(i,j)+vold(i,j-1)) * dx;
        
        % % combine UDS
        % ave = dve + max(-fve,0);
        % avw = dvw + max(fvw,0);
        % avn = dvn + max(-fvn,0);
        % avs = dvs + max(fvs,0);
        
        % avp = ave+avw+avn+avs + rho*dx*dy/dt;
        % AVp(i,j) = avp;
        
        % % RHS
        % rhs(idx) = -(pstar(i,j+1)-pstar(i,j))*dx + ...
            % rho*dx*dy*vold(i,j)/dt;
        
        % % Dirichlet BC
        % if (i == nx+1)
            % rhs(idx) = rhs(idx) + ave*vmac(i+1,j);
            % ave = 0;
        % end
        % if (i == 2)
            % rhs(idx) = rhs(idx) + avw*vmac(i-1,j);
            % avw = 0;
        % end
        % if (j == ny)
            % rhs(idx) = rhs(idx) + avn*vmac(i,j+1);
            % avn = 0;
        % end
        % if (j == 2)
            % rhs(idx) = rhs(idx) + avs*vmac(i,j-1);
            % avs = 0;
        % end
        
        % % relaxed equations
        % rhs(idx) = rhs(idx) + avp/relax*(1-relax)*vold(i,j);
        % if (j ~= 2)
            % H(idx,ids) = -avs;
        % end
        % if (i ~= 2)
            % H(idx,idw) = -avw;
        % end
        % H(idx,idx) = avp/relax;
        % if (i ~= nx+1)
            % H(idx,ide) = -ave;
        % end
        % if (j ~= ny)
            % H(idx,idn) = -avn;
        % end
    % end
% end

sol = H \ rhs;

% restore v*
vstar = vmac;
vstar(2:nx+1,2:ny) = reshape(sol,nx,ny-1);
% enforce BC
vstar = apply_vmac_bc(vstar);

return
end



