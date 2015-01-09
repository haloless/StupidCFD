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

% ## CIPAdvect2D

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-06

clc;
clear all;


% Zalesak disk
xlen = 1;
ylen = 1;

ncell = 100;
nx = ncell + 1;
ny = ncell + 1;
dx = xlen / (nx-1);
dy = ylen / (ny-1);

do_tan_tx = 1;
tan_coef = 0.9;


[X,Y] = ndgrid(linspace(0,xlen,nx),linspace(0,ylen,ny));
xc = 0.5;
yc = 0.5;
R = sqrt((X-xc).^2+(Y-yc).^2);
f0 = 1 * ((X-0.5).^2+(Y-0.75).^2<=0.17^2 & (abs(X-0.5)>0.03 | Y>0.85));
if (do_tan_tx)
    f0t = CIPTanTrans(f0,tan_coef);
else
    f0t = f0;
end
gx0 = zeros(nx,ny);
gy0 = zeros(nx,ny);
gx0(2:nx-1,2:ny-1) = 1/(2*dx) * (f0t(3:nx,2:ny-1) - f0t(1:nx-2,2:ny-1));
gy0(2:nx-1,2:ny-1) = 1/(2*dy) * (f0t(2:nx-1,3:ny) - f0t(2:nx-1,1:ny-2));

nstep = 800;
omega = 2*pi / nstep;
u = -omega * (Y-yc);
v = omega * (X-xc);

if (0)
    figure;
    % mesh(X',Y',f0');
    % surf(X',Y',f0');
    contourf(X',Y',f0',10); axis equal;
    hold on;
    quiver(X',Y',u',v');
    hold off;
    xlabel('X'); ylabel('Y');
    % figure;
end

f = f0;
gx = gx0;
gy = gy0;

I = 2:nx-1;
J = 2:ny-1;

dt = 1;
max_step = nstep;

for step = 1:max_step
    XX = -u(I,J) * dt;
    YY = -v(I,J) * dt;
    
    % us = sign(u(I,J)); us(find(us==0)) = 1;
    % vs = sign(v(I,J)); vs(find(vs==0)) = 1;
    % [Iup,Jup] = ndgrid(I,J);
    us = sign(u); us(find(us==0)) = 1;
    vs = sign(v); vs(find(vs==0)) = 1;
    [Iup,Jup] = ndgrid(1:nx,1:ny);
    Iup = Iup - us;
    Jup = Jup - vs;
    
    if (do_tan_tx)
        f = CIPTanTrans(f,tan_coef);
    end
    
    fn = f;
    gxn = gx;
    gyn = gy;
    
    A1 = zeros(nx,ny);
    E1 = zeros(nx,ny);
    B1 = zeros(nx,ny);
    F1 = zeros(nx,ny);
    D1 = zeros(nx,ny);
    C1 = zeros(nx,ny);
    G1 = zeros(nx,ny);
    
    for j = J
        for i = I
            iup = Iup(i,j);
            jup = Jup(i,j);
            zx = us(i,j);
            zy = vs(i,j);
            
            a1 = ((gx(iup,j)+gx(i,j))*dx*zx - 2*(f(i,j)-f(iup,j))) / (dx^3*zx);
            e1 = (3*(f(iup,j)-f(i,j)) + (gx(iup,j)+2*gx(i,j))*dx*zx) / dx^2;
            
            b1 = ((gy(i,jup)+gy(i,j))*dy*zy - 2*(f(i,j)-f(i,jup))) / (dy^3*zy);
            f1 = (3*(f(i,jup)-f(i,j)) + (gy(i,jup)+2*gy(i,j))*dy*zy) / dy^2;
            
            tmp = f(i,j) -f(i,jup) - f(iup,j) + f(iup,jup);
            tmq = gy(iup,j) - gy(i,j);
            
            d1 = (-tmp - tmq*dy*zy) / (dx*dy^2*zx);
            c1 = (-tmp - (gx(i,jup)-gx(i,j))*dx*zx) / (dx^2*dy*zy);
            g1 = (-tmq  + c1*dx^2) / (dx*zx);
            
            A1(i,j) = a1;
            E1(i,j) = e1;
            B1(i,j) = b1;
            F1(i,j) = f1;
            D1(i,j) = d1;
            C1(i,j) = c1;
            G1(i,j) = g1;
            
            % xx = XX(i,j);
            % yy = YY(i,j);
            % fn(i,j) = ((a1*xx+c1*yy+e1)*xx + g1*yy + gx(i,j)) * xx + ...
            % ((b1*yy+d1*xx+f1)*yy + gy(i,j)) * yy + f(i,j);
            % fx
        end
    end
    
    A1 = A1(I,J);
    E1 = E1(I,J);
    B1 = B1(I,J);
    F1 = F1(I,J);
    D1 = D1(I,J);
    C1 = C1(I,J);
    G1 = G1(I,J);
    
    fn(I,J) = ((A1.*XX+C1.*YY+E1).*XX + G1.*YY + gx(I,J)) .* XX + ...
        ((B1.*YY+D1.*XX+F1).*YY + gy(I,J)) .* YY + f(I,J);
    gxn(I,J) = (3*A1.*XX + 2*(C1.*YY+E1)).*XX + (D1.*YY+G1).*YY + gx(I,J);
    gyn(I,J) = (3*B1.*YY + 2*(D1.*XX+F1)).*YY + (C1.*XX+G1).*XX + gy(I,J);
    
    f = fn;
    gx = gxn;
    gy = gyn;
    
    if (do_tan_tx)
        f = CIPTanITrans(f,tan_coef);
    end
    
    if (mod(step,20)==0)
        subplot(1,2,1);
        contourf(X',Y',f',10); axis equal; 
        caxis([0 1]); colorbar;
        xlabel('X'); ylabel('Y');
        title(['contour ', int2str(step)]);
        
        subplot(1,2,2);
        mesh(X',Y',f');
        % surf(X',Y',f','EdgeColor','none');
        axis([0 xlen 0 ylen -0.1 1.1]);
        xlabel('X'); ylabel('Y');
        title(['disk ', int2str(step)]);
        
        drawnow;
    end
end





