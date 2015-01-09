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

% ## VelocityLapOp

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-26

function [ LapU,LapV ] = VelocityLapOp (nx,ny,dx,dy,dt)
%
EBGlobals;
%
coefx = nu*dt/2 * 1/dx^2;
coefy = nu*dt/2 * 1/dy^2;

% U
NU = (nx+1) * ny;
% LapU = spalloc(NU,NU,5*NU);
% % CorrU = spalloc(NU,NU,5*NU);
% corrU = zeros(NU,1);

% idx = 0;
% stride = nx + 1;
% for j = 1:ny
% for i = 1:nx+1
    % idx = idx + 1;
    % idw = idx - 1;
    % ide = idx + 1;
    % ids = idx - stride;
    % idn = idx + stride;
    
    % if (i==1)
        % LapU(idx,idx) = 1;
    % else
        % LapU(idx,idx) = 1;
        % % 
        % if (i~=2) % it is known that i~=1
            % LapU(idw,idx) = -coefx;
            % LapU(idx,idx) = LapU(idx,idx) + coefx;
        % else
            % LapU(idx,idx) = LapU(idx,idx) + coefx;
            % corrU(idx) = corrU(idx) + coefx*UIn;
        % end
        % %
        % if (j~=1)
            % LapU(ids,idx) = -coefy;
            % LapU(idx,idx) = LapU(idx,idx) + coefy;
        % else
            % LapU(idx,idx) = LapU(idx,idx) + 2*coefy;
        % end
        % %
        % if (i~=nx+1)
            % LapU(ide,idx) = -coefx;
            % LapU(idx,idx) = LapU(idx,idx) + coefx;
        % end
        % %
        % if (j~=ny)
            % LapU(idn,idx) = -coefy;
            % LapU(idx,idx) = LapU(idx,idx) + coefy;
        % else
            % LapU(idx,idx) = LapU(idx,idx) + 2*coefy;
        % end
    % end
% end
% end
% LapU = LapU';

% V
NV = nx * (ny+1);
% LapV = spalloc(NV,NV,5*NV);
% % CorrV = spalloc(NV,NV,5*NV);
% corrV = zeros(NV,1);

% idx = 0;
% stride = nx;
% for j = 1:ny+1
% for i = 1:nx
    % idx = idx + 1;
    % idw = idx - 1;
    % ide = idx + 1;
    % ids = idx - stride;
    % idn = idx + stride;
    
    % if (j==1 || j==ny+1)
        % LapV(idx,idx) = 1;
    % else
        % LapV(idx,idx) = 1;
        % % 1 < j < ny+1
        % if (j~=2)
            % LapV(ids,idx) = -coefy;
            % LapV(idx,idx) = LapV(idx,idx) + coefy;
        % else
            % LapV(idx,idx) = LapV(idx,idx) + coefy;
            % % corrV(idx) = corrV(idx) + coefy*0;
        % end
        % if (j~=ny)
            % LapV(idn,idx) = -coefy;
            % LapV(idx,idx) = LapV(idx,idx) + coefy;
        % else
            % LapV(idx,idx) = LapV(idx,idx) + coefy;
            % % corrV(idx) = corrV(idx) + coefy*0;
        % end
        % %
        % if (i~=1)
            % LapV(idw,idx) = -coefx;
            % LapV(idx,idx) = LapV(idx,idx) + coefx;
        % else
            % LapV(idx,idx) = LapV(idx,idx) + 2*coefx;
        % end
        % %
        % if (i~=nx)
            % LapV(ide,idx) = -coefx;
            % LapV(idx,idx) = LapV(idx,idx) + coefx;
        % end
    % end
% end
% end
% LapV = LapV';

% U
Aw = zeros(nx+1,ny);
Ae = zeros(nx+1,ny);
As = zeros(nx+1,ny);
An = zeros(nx+1,ny);
Ap = ones(nx+1,ny);
[I,J] = ndgrid(1:nx+1,1:ny);
Ind = I + (J-1)*(nx+1);
%
DirMask = (I==1); NonDirMask = ~DirMask;
% i-1
idx = (I~=1) & NonDirMask; idw = idx;
Aw(idx) = -coefx;
Aw(2,:) = 0;
Ap(idx) = Ap(idx) + coefx;
% j-1
idx = (J~=1) & NonDirMask; ids = idx;
As(idx) = -coefy;
Ap(idx) = Ap(idx) + coefy;
idx = (J==1) & NonDirMask;
Ap(idx) = Ap(idx) + coefy*2;
% j+1
idx = (J~=ny) & NonDirMask; idn = idx;
An(idx) = -coefy;
Ap(idx) = Ap(idx) + coefy;
idx = (J==ny) & NonDirMask;
Ap(idx) = Ap(idx) + coefy*2;
% i+1
idx = (I~=nx+1) & NonDirMask; ide = idx;
Ae(idx) = -coefx;
Ap(idx) = Ap(idx) + coefx;
idx = (I==nx+1) & NonDirMask;
Ap(idx) = Ap(idx) + 0;
%
tuples = [ ...
    reshape(Ind(idw),[],1), reshape(Ind(idw)-1,[],1), reshape(Aw(idw),[],1); ...
    reshape(Ind(ids),[],1), reshape(Ind(ids)-nx-1,[],1), reshape(As(ids),[],1); ...
    reshape(Ind(idn),[],1), reshape(Ind(idn)+nx+1,[],1), reshape(An(idn),[],1); ...
    reshape(Ind(ide),[],1), reshape(Ind(ide)+1,[],1), reshape(Ae(ide),[],1); ...
    reshape(Ind,[],1), reshape(Ind,[],1), reshape(Ap,[],1)];
LapU = sparse(tuples(:,1), tuples(:,2), tuples(:,3), NU,NU);

% V
Aw = zeros(nx,ny+1);
Ae = zeros(nx,ny+1);
As = zeros(nx,ny+1);
An = zeros(nx,ny+1);
Ap = ones(nx,ny+1);
[I,J] = ndgrid(1:nx,1:ny+1);
Ind = I + (J-1)*(nx);
%
DirMask = (J==1 | J==ny+1); NonDirMask = ~DirMask;
% i-1
idx = (I~=1) & NonDirMask; idw = idx;
Aw(idx) = -coefx;
Ap(idx) = Ap(idx) + coefx;
idx = (I==1) & NonDirMask;
Ap(idx) = Ap(idx) + coefx*2;
% j-1
idx = (J~=1) & NonDirMask; ids = idx;
As(idx) = -coefy;
As(:,2) = 0;
Ap(idx) = Ap(idx) + coefy;
% j+1
idx = (J~=ny+1) & NonDirMask; idn = idx;
An(idx) = -coefy;
An(:,ny) = 0;
Ap(idx) = Ap(idx) + coefy;
% i+1
idx = (I~=nx) & NonDirMask; ide = idx;
Ae(idx) = -coefx;
Ap(idx) = Ap(idx) + coefx;
idx = (I==nx) & NonDirMask;
Ap(idx) = Ap(idx) + 0;
%
tuples = [ ...
    reshape(Ind(idw),[],1), reshape(Ind(idw)-1,[],1), reshape(Aw(idw),[],1); ...
    reshape(Ind(ids),[],1), reshape(Ind(ids)-nx,[],1), reshape(As(ids),[],1); ...
    reshape(Ind(idn),[],1), reshape(Ind(idn)+nx,[],1), reshape(An(idn),[],1); ...
    reshape(Ind(ide),[],1), reshape(Ind(ide)+1,[],1), reshape(Ae(ide),[],1); ...
    reshape(Ind,[],1), reshape(Ind,[],1), reshape(Ap,[],1)];
LapV = sparse(tuples(:,1), tuples(:,2), tuples(:,3), NV,NV);


return
end
