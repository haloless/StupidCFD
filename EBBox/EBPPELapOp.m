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

% ## EBPPELapOp

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-29

function [ Lap,rhs_corr ] = EBPPELapOp (nx,ny,dx,dy, bx,by)
% Description

EBGlobals;

% check
if (size(bx,1)~=nx+1 || size(bx,2)~=ny)
    error('Beta_x wrong size');
end
if (size(by,1)~=nx || size(by,2)~=ny+1)
    error('Beta_y wrong size');
end

N = nx*ny;
Lap = spalloc(N,N,N*5);
rhs_corr = zeros(nx,ny);

cx = 1/dx^2;
cy = 1/dy^2;
tol = (cx+cy)/2 * 1e-4;
tol = 0;

I = 1:nx;
J = 1:ny;
Be = cx .* bx(I+1,J);
Bw = cx .* bx(I,J);
Bn = cy .* by(I,J+1);
Bs = cy .* by(I,J);

idx = 0;
for j = 1:ny
for i = 1:nx
    idx = idx + 1;
    idw = idx - 1;
    ide = idx + 1;
    idn = idx + nx;
    ids = idx - nx;
    
    Lkk = 0;
    if (i~=1 && abs(Bw(i,j))>tol)
        Lap(idw,idx) = -Bw(i,j);
        Lkk = Lkk + Bw(i,j);
        % Lap(idx,idx) = Lap(idx,idx) + cx*Bw(i,j);
    end
    if (j~=1 && abs(Bs(i,j))>tol)
        Lap(ids,idx) = -Bs(i,j);
        Lkk = Lkk + Bs(i,j);
        % Lap(idx,idx) = Lap(idx,idx) + cy*Bs(i,j);
    end
    if (j~=ny && abs(Bn(i,j))>tol)
        Lap(idn,idx) = -Bn(i,j);
        Lkk = Lkk + Bn(i,j);
        % Lap(idx,idx) = Lap(idx,idx) + cy*Bn(i,j);
    end
    if (i~=nx && abs(Be(i,j))>tol)
        Lap(ide,idx) = -Be(i,j);
        Lkk = Lkk + Be(i,j);
        % Lap(idx,idx) = Lap(idx,idx) + cx*Be(i,j);
    elseif (i==nx)
        Lkk = Lkk + Be(i,j);
        % Lap(idx,idx) = Lap(idx,idx) + cx*Be(i,j);
        rhs_corr(i,j) = rhs_corr(i,j) + Be(i,j)*POut;
    end
    
    if (abs(Lkk) <= tol*4)
        Lap(idx,idx) = 1;
    else
        Lap(idx,idx) = Lkk;
    end
end
end

Lap = Lap';
rhs_corr = reshape(rhs_corr,N,1);

return
end
