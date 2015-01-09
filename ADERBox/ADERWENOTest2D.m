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

% ## ADERWENOTest2D

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-01-19

clear all

% ader_order = 3;
% ader_order = 4;
ader_order = 5;

ADERWENOGlobals1D;
ADERWENOInit1D(ader_order);

M = MDegree;
N = NPoint;

eta = GausEta;

% test
xcen = 0.0;
ycen = 0.0;
hx = 0.01;
hy = 0.01;
% f = @(x,y) (x.^2 - y.^2);
f = @(x,y) (x.^2 + 3 .*x.*y - 2 .*y.^2);

nx = M*4 + 1;
ny = M*4 + 1;
icen = M + M + 1;
jcen = M + M + 1;
xs = linspace(xcen-2*M*hx, xcen+2*M*hx, nx);
ys = linspace(ycen-2*M*hy, xcen+2*M*hy, ny);
[xs,ys] = ndgrid(xs,ys);
if (0)
    fs = f(xs,ys);
else
    for j = jcen-M*2:jcen+M*2
    for i = icen-M*2:icen+M*2
        xlo = xs(i,j)-hx/2; xhi = xs(i,j)+hx/2;
        ylo = ys(i,j)-hy/2; yhi = ys(i,j)+hy/2;
        fs(i,j) = 1/(hx*hy) * integral2(f, xlo,xhi, ylo,yhi);
    end
    end
end

% reconstruction in x
qwx = zeros(nx,ny,N);
for j = jcen-M:jcen+M
for i = icen-M:icen+M
    qweno = ADERWENOReconstruct1D(fs(i-M:i+M,j), hx);
    qwx(i,j,:) = qweno;
end
end

% reconstruction in y
% central cell only
qwy = zeros(nx,ny,N,N);
for j = jcen:jcen
for i = icen:icen
    for l = 1:N
        ql = reshape(qwx(i,j-M:j+M,l), [],1);
        qweno = ADERWENOReconstruct1D(ql,hy);
        qwy(i,j,l,:) = qweno;
    end
end
end

qcs = reshape(qwy(icen,jcen,:,:), N, N);
xcs = (xs(icen,jcen)-hx/2) + eta*hx;
ycs = (ys(icen,jcen)-hy/2) * eta*hy;
[xcs,ycs] = ndgrid(xcs,ycs);

xes = (xs-(xcen-hx/2)) ./ hx;
yes = (ys-(ycen-hy/2)) ./ hy;
ves = zeros(nx,ny);
for j = 1:ny
for i = 1:nx
    for m = 1:N
    for l = 1:N
        vx = polyval(LagrPsi(:,l),xes(i,j));
        vy = polyval(LagrPsi(:,m),yes(i,j));
        ves(i,j) = ves(i,j) + qcs(l,m)*vx*vy;
    end
    end
end
end

if (1)
    figure;
    surf(xs,ys,fs);
    
    figure;
    mesh(xcs,ycs,qcs);
    hidden off;
    
    figure;
    surf(xs,ys,ves);
end



