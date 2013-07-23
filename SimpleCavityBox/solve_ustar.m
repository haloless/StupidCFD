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

% ## calc_u

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-02

function [ ustar ] = solve_ustar (umac,vmac,uold,vold,pstar)

simple_globals;

relax = urelax;

% ucell = zeros(nx+2,ny+2);
% ucell(2:nx+1,:) = 0.5 * (uold(1:nx,:) + uold(2:nx+1,:));

% vcell = zeros(nx+2,ny+2);
% vcell(:,2:ny+1) = 0.5 * (vold(:,1:ny) + vold(:,2:ny+1));


N = (nx-1) * ny;
A = spalloc(N,N,5*N);

Aw = zeros(nx+1,ny+2);
As = zeros(nx+1,ny+2);
An = zeros(nx+1,ny+2);
Ae = zeros(nx+1,ny+2);
Ap = zeros(nx+1,ny+2);

rhs = zeros(nx+1,ny+2);

% vectorized version
due = visc * dy / dx;
duw = visc * dy / dx;
dun = visc * dx / dy;
dus = visc * dx / dy;

fue = rho * dy * 0.5*(uold(3:nx+1,2:ny+1)+uold(2:nx,2:ny+1));
fuw = rho * dy * 0.5*(uold(1:nx-1,2:ny+1)+uold(2:nx,2:ny+1));
fun = rho * dx * 0.5*(vold(2:nx,2:ny+1)+vold(3:nx+1,2:ny+1));
fus = rho * dx * 0.5*(vold(2:nx,1:ny)+vold(3:nx+1,1:ny));

aue = zeros(nx+1,ny+2);
auw = zeros(nx+1,ny+2);
aun = zeros(nx+1,ny+2);
aus = zeros(nx+1,ny+2);

aue(2:nx,2:ny+1) = due + max(-fue,0);
auw(2:nx,2:ny+1) = duw + max(fuw,0);
aun(2:nx,2:ny+1) = dun + max(-fun,0);
aus(2:nx,2:ny+1) = dus + max(fus,0);
aup = aue+auw+aun+aus + rho*dx*dy/dt;

AUp(2:nx,2:ny+1) = aup(2:nx,2:ny+1);

rhs = -dy * diff(pstar) + rho*dx*dy/dt * uold;
% E
rhs(nx,2:ny+1) = rhs(nx,2:ny+1) + aue(nx,2:ny+1).*umac(nx+1,2:ny+1);
aue(nx,2:ny+1) = 0;
% W
rhs(2,2:ny+1) = rhs(2,2:ny+1) + auw(2,2:ny+1).*umac(1,2:ny+1);
auw(2,2:ny+1) = 0;
% N
rhs(2:nx,ny+1) = rhs(2:nx,ny+1) + aun(2:nx,ny+1).*umac(2:nx,ny+2);
aun(2:nx,ny+1) = 0;
% S
rhs(2:nx,2) = rhs(2:nx,2) + aus(2:nx,2).*umac(2:nx,1);
aus(2:nx,2) = 0;

rhs(2:nx,2:ny+1) = rhs(2:nx,2:ny+1) + (1-relax)/relax*aup(2:nx,2:ny+1).*uold(2:nx,2:ny+1);



idx = 0;
stride = nx-1;
for j = 2:ny+1
    for i = 2:nx
        idx = idx + 1;
        ide = idx + 1;
        idw = idx - 1;
        idn = idx + stride;
        ids = idx - stride;
        
        if (j ~= 2)
            A(idx,ids) = -aus(i,j);
        end
        if (i ~= 2)
            A(idx,idw) = -auw(i,j);
        end
        A(idx,idx) = aup(i,j) / relax;
        if (i ~= nx)
            A(idx,ide) = -aue(i,j);
        end
        if (j ~= ny+1)
            A(idx,idn) = -aun(i,j);
        end
        
    end
end

% idx = 0;
% stride = nx-1;
% for j = 2:ny+1
    % for i = 2:nx
        % idx = idx + 1;
        % ide = idx + 1;
        % idw = idx - 1;
        % idn = idx + stride;
        % ids = idx - stride;
        
        % % diffusion coef.
        % due = visc * dy / dx;
        % duw = visc * dy / dx;
        % dun = visc * dx / dy;
        % dus = visc * dx / dy;
        
        % % advection coef.
        % fue = rho * 0.5*(uold(i+1,j)+uold(i,j)) * dy;
        % fuw = rho * 0.5*(uold(i-1,j)+uold(i,j)) * dy;
        % fun = rho * 0.5*(vold(i,j)+vold(i+1,j)) * dx;
        % fus = rho * 0.5*(vold(i,j-1)+vold(i+1,j-1)) * dx;
        
        % aue = due + max(-fue,0);
        % auw = duw + max(fuw,0);
        % aun = dun + max(-fun,0);
        % aus = dus + max(fus,0);
        
        % % RHS
        % rhs(i,j) = -(pstar(i+1,j)-pstar(i,j))*dy + ...
            % rho*dx*dy*uold(i,j)/dt;
        
        % % BC
        % aup = aue+auw+aun+aus + rho*dx*dy/dt;
        % AUp(i,j) = aup;
        
        % if (i == nx)
            % rhs(i,j) = rhs(i,j) + aue*umac(i+1,j);
            % aue = 0;
        % end
        % if (i == 2)
            % rhs(i,j) = rhs(i,j) + auw*umac(i-1,j);
            % auw = 0;
        % end
        % if (j == ny+1)
            % rhs(i,j) = rhs(i,j) + aun*umac(i,j+1);
            % aun = 0;
        % end
        % if (j == 2)
            % rhs(i,j) = rhs(i,j) + aus*umac(i,j-1);
            % aus = 0;
        % end
        
        % % relaxed equations
        % Ap(i,j) = aup / urelax;
        % Ae(i,j) = -aue;
        % Aw(i,j) = -auw;
        % An(i,j) = -aun;
        % As(i,j) = -aus;
        % rhs(i,j) = rhs(i,j) + aup/urelax*(1-urelax)*uold(i,j);
        
        % if (j ~= 2)
            % A(idx,ids) = As(i,j);
        % end
        % if (i ~= 2)
            % A(idx,idw) = Aw(i,j);
        % end
        % A(idx,idx) = Ap(i,j);
        % if (i ~= nx)
            % A(idx,ide) = Ae(i,j);
        % end
        % if (j ~= ny+1)
            % A(idx,idn) = An(i,j);
        % end
    % end
% end

rhs = reshape(rhs(2:nx,2:ny+1),N);

% solve
sol = A \ rhs;

% get u*
ustar = umac;
ustar(2:nx,2:ny+1) = reshape(sol,nx-1,ny);
% enforce BC
ustar = apply_umac_bc(ustar);

return
end



