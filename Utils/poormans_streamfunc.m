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

% ## poormans_streamfunc

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-03

function [ psi ] = poormans_streamfunc (nx,ny,dx,dy,umac,vmac)

% Description:
% Note that Psi is defined on cell nodes!
% u = d(psi)/dy
% v = -d(psi)/dx
% a Poisson equation is solved 
% with Dirichlet BC (psi=0 on boundary)


Npsi = (nx-1) * (ny-1);

Lpsi = kron(speye(ny-1),K1(nx-1,dx,2))+kron(K1(ny-1,dy,2),speye(nx-1));

rhs = zeros(nx-1,ny-1);
rhs = rhs - 1/dy * (umac(2:nx,3:ny+1) - umac(2:nx,2:ny));
rhs = rhs + 1/dx * (vmac(3:nx+1,2:ny) - vmac(2:nx,2:ny));
rhs = reshape(rhs,Npsi);

perm = symamd(Lpsi);
RLpsi = chol(Lpsi(perm,perm));
% RLpsit = RLpsi';

psi = zeros(1,Npsi);
psi(perm) = RLpsi \ (RLpsi' \ rhs(perm));

% psi = Lpsi \ rhs;
psi = reshape(psi,nx-1,ny-1);

return
end
