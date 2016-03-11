% ## Copyright (C) 2014 homu
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

% ## TestMG1D

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-01-25


function [ Ac Aw Ae As An rhs ] = TestInitProb(nx,ny,dx,dy)

%
gravy = -1;

% bc(dir,side)
% probbc = [ 1.0, 1.0; 1.0, 1.0 ];
probbc = [ 1.0, 1.0; 1.0, -1.0 ];

% setup test problem
xcell = zeros(nx,ny);
ycell = zeros(nx,ny);
for j = 1:ny
for i = 1:nx
    xcell(i,j) = (i - 0.5) * dx;
    ycell(i,j) = (j - 0.5) * dy;
end
end

%
dens = zeros(nx,ny);
dens(:) = 1.0;
dens(ycell(:)<=0.5) = 1.0e3;

% build rhs
rhs = zeros(nx,ny);
for j = 1:ny
for i = 1:nx
    vh = gravy;
    if (j == ny) 
        vh = 0; 
    end
    
    vl = gravy;
    if (j == 1) 
        vl = 0; 
    end
    
    divu = (vh-vl)/dy;
    rhs(i,j) = -divu;
end
end

% build matrix
Ac = zeros(nx,ny);
Aw = zeros(nx,ny);
Ae = zeros(nx,ny);
As = zeros(nx,ny);
An = zeros(nx,ny);
for j = 1:ny
for i = 1:nx
    rhow = dens(i,j);
    if (i > 1)
        rhow = 0.5 * (dens(i,j)+dens(i-1,j));
    end
    rhoe = dens(i,j);
    if (i < nx)
        rhoe = 0.5 * (dens(i,j)+dens(i+1,j));
    end
    rhos = dens(i,j);
    if (j > 1)
        rhos = 0.5 * (dens(i,j)+dens(i,j-1));
    end
    rhon = dens(i,j);
    if (j < ny)
        rhon = 0.5 * (dens(i,j)+dens(i,j+1));
    end
    
    %
    axlo = 1.0/rhow/(dx*dx);
    axhi = 1.0/rhoe/(dx*dx);
    aylo = 1.0/rhos/(dy*dy);
    ayhi = 1.0/rhon/(dy*dy);
    
    %
    bxlo = 0;
    if (i > 1)
        bxlo = -axlo;
    end
    bxhi = 0;
    if (i < nx)
        bxhi = -axhi;
    end
    bylo = 0;
    if (j > 1)
        bylo = -aylo;
    end
    byhi = 0;
    if (j < ny)
        byhi = -ayhi;
    end
    
    %
    cxlo = probbc(1,1) * axlo;
    if (i > 1)
        cxlo = 0;
    end
    cxhi = probbc(1,2) * axhi;
    if (i < nx)
        cxhi = 0;
    end
    cylo = probbc(2,1) * aylo;
    if (j > 1)
        cylo = 0;
    end
    cyhi = probbc(2,2) * ayhi;
    if (j < ny)
        cyhi = 0;
    end
    
    %
    Aw(i,j) = bxlo;
    Ae(i,j) = bxhi;
    As(i,j) = bylo;
    An(i,j) = byhi;
    Ac(i,j) = (axlo+axhi+aylo+ayhi) - (cxlo+cxhi+cylo+cyhi);
end
end


return
end






