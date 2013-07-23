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

% ## WENO_LuOp

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-16

function [ Lu ] = WENO_LuOp (u, F, dFdu, nx, ngrow, dx)
% Description:

% WENO_BurgersGlobals;

    % apply BC
    u = WENO_bc(u, nx, ngrow);
    
    % flux splitting
    [fp,fn] = WENO_fluxsplit(u,F,dFdu);
    
    % reconstruct flux at cells interface
    hp = zeros(size(fp));
    hn = zeros(size(fn));
    for i = 1+ngrow:nx-ngrow
        xr = i-ngrow:i+ngrow;
        % [h_right,h_left] = WENO_flux(fp(xr),fn(xr));
        % hn(i) = h_right;
        % hp(i-1) = h_left;
        [hn(i),hp(i-1)] = WENO_flux(fp(xr),fn(xr));
    end
    
    % f_1/2 = f_1/2^plus + f_1/2^minus
    h = hp + hn;
    
    Lu = zeros(size(h));
    Lu(2:end) = -1/dx * diff(h);
    % h_right = h;
    % h_left = [0, h(1:end-1)];
    
    % update
    % u_new = u - dt/dx * (h_right-h_left);
    % u_new(2:end) = u(2:end) - dt/dx * diff(h);
    % u_new = WENO_bc(u_new, nx, ngrow);
return
end


