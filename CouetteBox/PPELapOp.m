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

% ## PPELapOp

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-08

function [ Lap,rhs_corr ] = PPELapOp (nx,ny,dx,dy)

EBGlobals;

N = nx*ny;

Lap = sparse(N,N);
rhs_corr = zeros(nx,ny);


if (0)
    warning('Obsolete routine');
    
    Lap = sparse(nx*ny,nx*ny);
    idx = 0;
    for j = 1:ny
        for i = 1:nx
            idx = idx + 1;
            idw = idx - 1;
            ide = idx + 1;
            ids = idx - nx;
            idn = idx + nx;
            
            if (i ~= 1)
                Lap(idx,idw) = -1/dx^2;
                Lap(idx,idx) = Lap(idx,idx) + 1/dx^2;
            end
            if (j ~= 1)
                Lap(idx,ids) = -1/dy^2;
                Lap(idx,idx) = Lap(idx,idx) + 1/dy^2;
            end
            if (j ~= ny)
                Lap(idx,idn) = -1/dy^2;
                Lap(idx,idx) = Lap(idx,idx) + 1/dy^2;
            end
            if (i ~= nx)
                Lap(idx,ide) = -1/dx^2;
                Lap(idx,idx) = Lap(idx,idx) + 1/dx^2;
            else
                Lap(idx,idx) = Lap(idx,idx) + 1/dx^2;
            end
        end
    end
else
    Aw = zeros(nx,ny);
    Ae = zeros(nx,ny);
    As = zeros(nx,ny);
    An = zeros(nx,ny);
    Ap = zeros(nx,ny);
    [Ip,Jp] = ndgrid(1:nx,1:ny);
    Id = Ip + (Jp-1)*nx;


    idx = (Ip ~= 1);
    idw = idx;
    Aw(idx) = -1/dx^2;
    Ap(idx) = Ap(idx) + 1/dx^2;

    idx = (Jp ~= 1);
    ids = idx;
    As(idx) = -1/dy^2;
    Ap(idx) = Ap(idx) + 1/dy^2;

    idx = (Jp ~= ny);
    idn = idx;
    An(idx) = -1/dy^2;
    Ap(idx) = Ap(idx) + 1/dy^2;

    idx = (Ip ~= nx);
    ide = idx;
    Ae(idx) = -1/dx^2;
    Ap(idx) = Ap(idx) + 1/dx^2;
    if (0)
    % Dirichlet
    ridx = ~idx;
    Ap(ridx) = Ap(ridx) + 1/dx^2;
    rhs_corr(ridx) = rhs_corr(ridx) + 1/dx^2 * POut;
    end

    tuples = [ ...
        reshape(Id(idw),[],1), reshape(Id(idw)-1,[],1), reshape(Aw(idw),[],1); ...
        reshape(Id(ids),[],1), reshape(Id(ids)-nx,[],1), reshape(As(ids),[],1); ...
        reshape(Id(idn),[],1), reshape(Id(idn)+nx,[],1), reshape(An(idn),[],1); ...
        reshape(Id(ide),[],1), reshape(Id(ide)+1,[],1), reshape(Ae(ide),[],1); ...
        reshape(Id,[],1),      reshape(Id,[],1),      reshape(Ap,[],1)];
    Lap = sparse(tuples(:,1),tuples(:,2),tuples(:,3),N,N);

end

rhs_corr = reshape(rhs_corr,[],1);

return
end
