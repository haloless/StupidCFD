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

function [ Lap ] = PMVelocityLapOp (eb_vof,nx,ny,dx,dy,dt)
%
EBGlobals;

% compute (penalty) viscosity from harmonic average
[nucell,nunode] = PMViscosity(eb_vof,nx,ny);

%
dx2 = dx^2;
dy2 = dy^2;
dxdy = dx * dy;
%
% uind = zeros(nx+3,ny+2);
% uind(3:nx+1,2:ny+1) = reshape((1:(nx-1)*ny), nx-1,ny);
% vind = zeros(nx+2,ny+3);
% vind(2:nx+1,3:ny+1) = reshape((1:nx*(ny-1))+(nx-1)*ny, nx,ny-1);
uind = zeros(nx+1,ny+2);
uind(2:nx,2:ny+1) = reshape((1:(nx-1)*ny), nx-1,ny);
vind = zeros(nx+2,ny+1);
vind(2:nx+1,2:ny) = reshape((1:nx*(ny-1))+(nx-1)*ny, nx,ny-1);


%
N = (nx-1)*ny + nx*(ny-1);
nzmax = 9 * N;
is = zeros(nzmax,1);
js = zeros(nzmax,1);
ss = zeros(nzmax,1);
% sr = zeros(nzmax,1); % Dirichlet BC contribution

%
idx = 1;

%
for j = 2:ny+1
for i = 2:nx
    % diagonal
    is(idx) = uind(i,j);
    js(idx) = uind(i,j);
    ss(idx) = 1.0;
    
    %
    coef = dt * 2.0 * nucell(i,j) / dx2;
    ss(idx) = ss(idx) + coef;
    if uind(i-1,j) > 0
        is(idx+1) = uind(i,j);
        js(idx+1) = uind(i-1,j);
        ss(idx+1) = -coef;
    else
        is(idx+1) = uind(i,j);
        js(idx+1) = uind(i,j);
        ss(idx+1) = 0.0;
    end
    %
    coef = dt * 2.0 * nucell(i+1,j) / dx2;
    ss(idx) = ss(idx) + coef;
    if uind(i+1,j) > 0
        is(idx+2) = uind(i,j);
        js(idx+2) = uind(i+1,j);
        ss(idx+2) = -coef;
    else
        is(idx+2) = uind(i,j);
        js(idx+2) = uind(i,j);
        ss(idx+2) = 0.0;
    end
    %
    coef = dt * nunode(i,j-1) / dy2;
    ss(idx) = ss(idx) + coef;
    if uind(i,j-1) > 0
        is(idx+3) = uind(i,j);
        js(idx+3) = uind(i,j-1);
        ss(idx+3) = -coef;
    else
        is(idx+3) = uind(i,j);
        js(idx+3) = uind(i,j);
        ss(idx+3) = coef;
    end
    %
    coef = dt * nunode(i,j) / dy2;
    ss(idx) = ss(idx) + coef;
    if uind(i,j+1) > 0
        is(idx+4) = uind(i,j);
        js(idx+4) = uind(i,j+1);
        ss(idx+4) = -coef;
    else
        is(idx+4) = uind(i,j);
        js(idx+4) = uind(i,j);
        ss(idx+4) = coef;
    end
    
    %
    coef = dt * nunode(i,j-1) / dxdy;
    if vind(i,j-1) > 0
        is(idx+5) = uind(i,j);
        js(idx+5) = vind(i,j-1);
        ss(idx+5) = -coef;
    else
        is(idx+5) = uind(i,j);
        js(idx+5) = uind(i,j);
        ss(idx+5) = 0.0;
    end
    if vind(i+1,j-1) > 0
        is(idx+6) = uind(i,j);
        js(idx+6) = vind(i+1,j-1);
        ss(idx+6) = coef;
    else
        is(idx+6) = uind(i,j);
        js(idx+6) = uind(i,j);
        ss(idx+6) = 0.0;
    end
    %
    coef = dt * nunode(i,j) / dxdy;
    if vind(i,j) > 0
        is(idx+7) = uind(i,j);
        js(idx+7) = vind(i,j);
        ss(idx+7) = coef;
    else
        is(idx+7) = uind(i,j);
        js(idx+7) = uind(i,j);
        ss(idx+7) = 0.0;
    end
    if vind(i+1,j) > 0
        is(idx+8) = uind(i,j);
        js(idx+8) = vind(i+1,j);
        ss(idx+8) = -coef;
    else
        is(idx+8) = uind(i,j);
        js(idx+8) = uind(i,j);
        ss(idx+8) = 0.0;
    end
    
    %
    idx = idx + 9;
end
end

%
for j = 2:ny
for i = 2:nx+1
    % diagonal
    is(idx) = vind(i,j);
    js(idx) = vind(i,j);
    ss(idx) = 1.0;
    
    %
    coef = dt * nunode(i-1,j) / dx2;
    ss(idx) = ss(idx) + coef;
    if vind(i-1,j) > 0
        is(idx+1) = vind(i,j);
        js(idx+1) = vind(i-1,j);
        ss(idx+1) = -coef;
    else
        is(idx+1) = vind(i,j);
        js(idx+1) = vind(i,j);
        ss(idx+1) = coef;
    end
    %
    coef = dt * nunode(i,j) / dx2;
    ss(idx) = ss(idx) + coef;
    if vind(i+1,j) > 0
        is(idx+2) = vind(i,j);
        js(idx+2) = vind(i+1,j);
        ss(idx+2) = -coef;
    else
        is(idx+2) = vind(i,j);
        js(idx+2) = vind(i,j);
        ss(idx+2) = coef;
    end
    
    %
    coef = dt * 2.0 * nucell(i,j) / dy2;
    ss(idx) = ss(idx) + coef;
    if vind(i,j-1) > 0
        is(idx+3) = vind(i,j);
        js(idx+3) = vind(i,j-1);
        ss(idx+3) = -coef;
    else
        is(idx+3) = vind(i,j);
        js(idx+3) = vind(i,j);
        ss(idx+3) = 0.0;
    end
    %
    coef = dt * 2.0 * nucell(i,j+1) / dy2;
    ss(idx) = ss(idx) + coef;
    if vind(i,j+1) > 0
        is(idx+4) = vind(i,j);
        js(idx+4) = vind(i,j+1);
        ss(idx+4) = -coef;
    else
        is(idx+4) = vind(i,j);
        js(idx+4) = vind(i,j);
        ss(idx+4) = 0.0;
    end
    
    %
    coef = dt * nunode(i-1,j) / dxdy;
    if uind(i-1,j) > 0
        is(idx+5) = vind(i,j);
        js(idx+5) = uind(i-1,j);
        ss(idx+5) = -coef;
    else
        is(idx+5) = vind(i,j);
        js(idx+5) = vind(i,j);
        ss(idx+5) = 0.0;
    end
    if uind(i-1,j+1) > 0
        is(idx+6) = vind(i,j);
        js(idx+6) = uind(i-1,j+1);
        ss(idx+6) = coef;
    else
        is(idx+6) = vind(i,j);
        js(idx+6) = vind(i,j);
        ss(idx+6) = 0.0;
    end
    %
    coef = dt * nunode(i,j) / dxdy;
    if uind(i,j) > 0
        is(idx+7) = vind(i,j);
        js(idx+7) = uind(i,j);
        ss(idx+7) = coef;
    else
        is(idx+7) = vind(i,j);
        js(idx+7) = vind(i,j);
        ss(idx+7) = 0.0;
    end
    if uind(i,j+1) > 0
        is(idx+8) = vind(i,j);
        js(idx+8) = uind(i,j+1);
        ss(idx+8) = -coef;
    else
        is(idx+8) = vind(i,j);
        js(idx+8) = vind(i,j);
        ss(idx+8) = 0.0;
    end
    
    %
    idx = idx + 9;
end
end

% find(is==0)
% find(js==0)

%
Lap = sparse(is,js,ss,N,N,nzmax);


return
end
