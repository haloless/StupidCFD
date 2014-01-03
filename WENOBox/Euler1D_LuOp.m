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

% ## Euler1D_LuOp

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-01-03

function [ Lu,dtau ] = Euler1D_LuOp (ucons, lo,hi,nx,ng,dx)

Euler1D_globals;

% conservative -> primitive
uprim = zeros(size(ucons));
uprim = Euler1D_ConsToPrim(ucons,uprim, lo,hi);

% apply BC
for i = 1:ng
    ucons(:,lo-i) = ucons(:,lo+i-1);
    ucons(:,hi+i) = ucons(:,hi-i+1);
    uprim(:,lo-i) = uprim(:,lo+i-1);
    uprim(:,hi+i) = uprim(:,hi-i+1);
end


% cell flux
fcell = Euler1D_flux(ucons, uprim);
% sound speed
sound = Euler1D_csound(uprim);
sig = zeros(size(sound));

% numerical flux
fhat = zeros(3,nx+1);

for i = lo:hi+1 % loop face
    % local LF flux splitting
    js =  i-3:i+2;
    alpha = max(abs(uprim(QVX,js)) + sound(js));
    fps = 0.5 .* (fcell(:,js) + alpha.*ucons(:,js));
    fms = 0.5 .* (fcell(:,js) - alpha.*ucons(:,js));
    % save local speed
    sig(i) = alpha;
    
    % Roe's mean matrix
    [Rhalf, Dhalf, Lhalf] = Euler1D_RoeMatrix(ucons(:,i-1),ucons(:,i),uprim(:,i-1),uprim(:,i));
    % transform to local characteristic fields
    ups = Lhalf * fps;
    ums = Lhalf * fms;
    
    % WENO resconstruction on each component
    upc = zeros(3,1);
    umc = zeros(3,1);
    for comp = 1:3
        upc(comp) = WENO3_upwind(ups(comp,:)');
        umc(comp) = WENO3_upwind(flipud(ums(comp,:)'));
    end
    
    % transform back
    fpc = Rhalf * upc;
    fmc = Rhalf * umc;
    
    % form the flux
    fhat(:,i) = fpc + fmc;
end

% conservative flux
Lu = zeros(size(ucons));
irange = lo:hi;
Lu(:,irange) = 1/dx * (fhat(:,irange)-fhat(:,irange+1));

% max time step allowable
dtau = dx / (max(sig) + 1.0e-8);

return
end



