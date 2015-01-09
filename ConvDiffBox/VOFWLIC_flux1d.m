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

% ## VOFWLIC_flux1d

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-10-04

function [ fs ] = VOFWLIC_flux1d (ncell,fs,us,ns,dt,dh,dir,ndim)
% Description
% ns: ncell*ndim

eps_small = 1e-10;

% THINC transformation
beta = 3.5;

% flux at cell face
flux = zeros(ncell+3,1);

nabs = abs(ns);
nsum = sum(nabs,2);

usign = (us > 0);

% interfacial cells
fflag = ((fs>eps_small) & (fs<1.0-eps_small));
nflag = (nsum > eps_small);

% looping edges
for i = 2:ncell+2
    if (us(i) <= 0)
        ii = i;
    else
        ii = i-1;
    end
    nn = nsum(ii);
    
    uu = us(i) * dt;
    if (abs(uu) > dh)
        error('dt too large, CFL condition hit');
    end
    
    if (~fflag(ii) || ~nflag(ii))
        % full liquid/gas cell, upwind
        flux(i) = fs(ii) * uu;
    else
        ilo = max(ii-1,1);
        ihi = min(ii+1,ncell+2);
        % direction of interface
        if (fs(ilo)<=fs(ihi))
            alpha = 1;
        else
            alpha = -1;
        end
        % THINC 1D
        a1 = exp(beta/alpha * (2*fs(ii)-1));
        a3 = exp(beta);
        xc = 0.5/beta * log(a3*(a3-a1)/(a1*a3-1));
        a4 = cosh(beta * (usign(i) - uu/dh - xc));
        a5 = cosh(beta * (usign(i) - xc));
        %
        F = 0.5 * (uu - alpha*dh/beta * log(a4/a5));
        w = nabs(ii,dir) / nn;
        flux(i) = w*F + (1-w)*fs(ii)*uu;
    end
end

% BC
% flux(1) = 0;
% flux(2) = 0;
% flux(ncell+1) = 0;
% flux(ncell+2) = 0;

I = 2:ncell+1;
fs(I) = fs(I) - 1/dh * (flux(I+1)-flux(I));

% BC
fs(1) = fs(2);
fs(ncell+2) = fs(ncell+1);

return
end
