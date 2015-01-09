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

% ## ../LevelsetBox/LSTest2

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-15


function LSTest2

clc;
clear all;

nx = 128 + 1;
ny = 128 + 1;

dx = 1 / (nx-1);
dy = 1 / (ny-1);
dh = min([dx dy]);

dt = 4e-3;
dtau = dt*1;

nreinit = 1000;


xs = linspace(0,1,nx);
ys = linspace(0,1,ny);


[X,Y] = ndgrid(xs,ys);

% initial condition
phi = sqrt((X-0.45).^2+((Y-0.8)*3).^2) - 0.4;

if (1)
    gx0 = zeros(nx,ny);
    gy0 = zeros(nx,ny);
    % gx0(:) = NaN;
    % gy0(:) = NaN;
    I = 2:nx-1; J = 2:ny-1;
    gx0(I,J) = 1/(2*dx) * (phi(I+1,J) - phi(I-1,J));
    gy0(I,J) = 1/(2*dy) * (phi(I,J+1) - phi(I,J-1));
    
    figure;
    subplot(2,2,1);
    surf(X',Y',phi', 'EdgeAlpha',0.1);
    subplot(2,2,2);
    surf(X',Y',gx0', 'EdgeAlpha',0.1);
    subplot(2,2,3);
    surf(X',Y',gy0', 'EdgeAlpha',0.1);
    
    figure;
    
    subplot(1,2,1)
    % surf(X',Y',phi','EdgeAlpha',0.2);
    [cont,ch] = contour(X',Y',phi');
    clabel(cont,ch);
    axis equal;
    axis([0 1 0 1]);
    title('Raw \phi')
    
    subplot(1,2,2);

    surf(X',Y',sqrt(gx0.^2+gy0.^2)', 'EdgeAlpha',0.1);
    title('\nabla \phi');
end

for ir = 1:nreinit
    phi_corr = dtau*FabsgradP(phi, dh, phi./sqrt(phi.^2+(2*dh)^2), 1);
    phi = phi - phi_corr;
    resid = norm(phi_corr(:),1);
end

if (1)
    figure;
    subplot(1,2,1);
    % figure;
    % surf(X',Y',phi','EdgeAlpha',0.2);
    [cont,ch] = contour(X',Y',phi');
    clabel(cont,ch);
    axis equal;
    axis([0 1 0 1]);
    title(['reinitializing:', int2str(ir),'/',int2str(nreinit), ...
    ';residual=', num2str(resid)]);
    
    subplot(1,2,2);
    gx0 = zeros(nx,ny);
    gy0 = zeros(nx,ny);
    gx0(:) = NaN;
    gy0(:) = NaN;
    I = 2:nx-1; J = 2:ny-1;
    gx0(I,J) = 1/(2*dx) * (phi(I+1,J) - phi(I-1,J));
    gy0(I,J) = 1/(2*dy) * (phi(I,J+1) - phi(I,J-1));
    surf(X',Y',sqrt(gx0.^2+gy0.^2)', 'EdgeAlpha',0.1);
    title('\nabla \phi');
   
    
end

end % end of LSTest2

%===================

function dP = FabsgradP(P,h,F,c)
if nargin<4, c = 0; if nargin<3, F = 1; end, end
DxP = diff(P)/h;   DxmP = DxP([1 1:end],:); DxpP = DxP([1:end end],:);
DyP = diff(P')'/h; DymP = DyP(:,[1 1:end]); DypP = DyP(:,[1:end end]);
Np = sqrt(max(DxmP,0).^2+min(DxpP,0).^2+max(DymP,0).^2+min(DypP,0).^2);
Nm = sqrt(min(DxmP,0).^2+max(DxpP,0).^2+min(DymP,0).^2+max(DypP,0).^2);
dP = max(F,0).*(Np-c)+min(F,0).*(Nm-c);
end

