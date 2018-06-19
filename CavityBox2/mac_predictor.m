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

% ## mac_residual_corr_form

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-06-26

function [ Hx,Hy ] = mac_predictor (umac,vmac, nx,ny, dx,dy, dt, Re, BC)
% Description:

% blend factor of CDS and UDS for advection
% gamma=0: fully CDS
% gamma=1: fully UDS
% umax = norm(umac(:),Inf);
% vmax = norm(vmac(:),Inf);
% gamma = 1.2 * max(umax*dt/dx, vmax*dt/dy);
% gamma = max(min(gamma,1),0);

dx2 = dx^2;
dy2 = dy^2;

Hx = zeros(nx+2,ny+2);

xrange = 2:nx+1;
yrange = 2:ny+1;

ue = 0.5 * (umac(xrange+1,yrange) + umac(xrange,yrange));
uw = 0.5 * (umac(xrange-1,yrange) + umac(xrange,yrange));
un = 0.5 * (umac(xrange,yrange+1) + umac(xrange,yrange));
us = 0.5 * (umac(xrange,yrange-1) + umac(xrange,yrange));
vn = 0.5 * (vmac(xrange-1,yrange+1) + vmac(xrange,yrange+1));
vs = 0.5 * (vmac(xrange-1,yrange) + vmac(xrange,yrange));

due = 0.5 * (umac(xrange+1,yrange) - umac(xrange,yrange));
duw = 0.5 * (umac(xrange,yrange) - umac(xrange-1,yrange));
dun = 0.5 * (umac(xrange,yrange+1) - umac(xrange,yrange));
dus = 0.5 * (umac(xrange,yrange) - umac(xrange,yrange-1));

% % hmm... not bad
uadv_cds = 1/dx*(ue.^2-uw.^2) + 1/dy*(un.*vn-us.*vs);
% uadv_uds = 1/dx * max(uw,0) .* (umac(xrange,yrange)-umac(xrange-1,yrange)) + ...
    % 1/dx * max(-ue,0) .* (umac(xrange,yrange)-umac(xrange+1,yrange)) + ...
    % 1/dy * max(vs,0) .* (umac(xrange,yrange)-umac(xrange,yrange-1)) + ...
    % 1/dy * max(-vn,0) .* (umac(xrange,yrange)-umac(xrange,yrange+1));
% uadv = (1-gamma)*uadv_cds + gamma*uadv_uds;

% % mixed CDS and UDS
% uadv = 1/dx * (ue.^2 - gamma*abs(ue).*due) - ...
    % 1/dx * (uw.^2 - gamma*abs(uw).*duw) + ...
    % 1/dy * (un.*vn - gamma*abs(vn).*dun) - ...
    % 1/dy * (us.*vs - gamma*abs(vs).*dus);

% 2nd-order UDS
% uadv = 1/dx * (max(ue,0).*umac(xrange,yrange) + min(ue,0).*umac(xrange+1,yrange)) - ...
    % 1/dx * (max(uw,0).*umac(xrange-1,yrange) + min(uw,0).*umac(xrange,yrange)) + ...
    % 1/dy * (max(vn,0).*umac(xrange,yrange) + min(vn,0).*umac(xrange,yrange+1)) - ...
    % 1/dy * (max(vs,0).*umac(xrange,yrange-1) + min(vs,0).*umac(xrange,yrange));
uadv = uadv_cds;
% uadv = 0;

d2udx2 = 1/dx2 * (umac(xrange+1,yrange) - 2*umac(xrange,yrange) + umac(xrange-1,yrange));
d2udy2 = 1/dy2 * (umac(xrange,yrange+1) - 2*umac(xrange,yrange) + umac(xrange,yrange-1));
Hx(xrange,yrange) = -uadv + 1/Re*(d2udx2+d2udy2);


Hy = zeros(nx+2,ny+2);

xrange = 2:nx+1;
yrange = 2:ny+1;

ue = 0.5 * (umac(xrange+1,yrange-1) + umac(xrange+1,yrange));
ve = 0.5 * (vmac(xrange+1,yrange) + vmac(xrange,yrange));
uw = 0.5 * (umac(xrange,yrange-1) + umac(xrange,yrange));
vw = 0.5 * (vmac(xrange-1,yrange) + vmac(xrange,yrange));
vn = 0.5 * (vmac(xrange,yrange+1) + vmac(xrange,yrange));
vs = 0.5 * (vmac(xrange,yrange-1) + vmac(xrange,yrange));

dve = 0.5 * (vmac(xrange+1,yrange) - vmac(xrange,yrange));
dvw = 0.5 * (vmac(xrange,yrange) - vmac(xrange-1,yrange));
dvn = 0.5 * (vmac(xrange,yrange+1) - vmac(xrange,yrange));
dvs = 0.5 * (vmac(xrange,yrange) - vmac(xrange,yrange-1));

vadv_cds = 1/dx*(ue.*ve-uw.*vw) + 1/dy*(vn.^2-vs.^2);
% vadv_uds = 1/dx * max(uw,0) .* (vmac(xrange,yrange)-vmac(xrange-1,yrange)) + ...
    % 1/dx * max(-ue,0) .* (vmac(xrange,yrange)-vmac(xrange+1,yrange)) + ...
    % 1/dy * max(vs,0) .* (vmac(xrange,yrange)-vmac(xrange,yrange-1)) + ...
    % 1/dy * max(-vn,0) .* (vmac(xrange,yrange)-vmac(xrange,yrange+1));
% vadv = (1-gamma)*vadv_cds + gamma*vadv_uds;

% % mixed CDS and UDS
% vadv = 1/dx * (ue.*ve - gamma*abs(ue).*dve) - ...
    % 1/dx * (uw.*vw - gamma*abs(uw).*dvw) + ...
    % 1/dy * (vn.^2 - gamma*abs(vn).*dvn) - ...
    % 1/dy * (vs.^2 - gamma*abs(vs).*dvs);

% 2nd-order UDS
% vadv = 1/dx * (max(ue,0).*vmac(xrange,yrange) + min(ue,0).*vmac(xrange+1,yrange)) - ...
    % 1/dx * (max(uw,0).*vmac(xrange-1,yrange) + min(uw,0).*vmac(xrange,yrange)) + ...
    % 1/dy * (max(vn,0).*vmac(xrange,yrange) + min(vn,0).*vmac(xrange,yrange+1)) - ...
    % 1/dy * (max(vs,0).*vmac(xrange,yrange-1) + min(vs,0).*vmac(xrange,yrange));
vadv = vadv_cds;
% vadv = 0;

d2vdx2 = 1/dx2 * (vmac(xrange+1,yrange) - 2*vmac(xrange,yrange) + vmac(xrange-1,yrange));
d2vdy2 = 1/dy2 * (vmac(xrange,yrange+1) - 2*vmac(xrange,yrange) + vmac(xrange,yrange-1));
Hy(xrange,yrange) = -vadv + 1/Re*(d2vdx2+d2vdy2);

return
end


