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

% ## Vorticity_incr

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-18

function [ dw_hat ] = Vorticity_incr (w_hat, kx,ky, nu,k2_visc,k2_ppe, cutoff_filter)

% stream function
psi_hat = -w_hat ./ k2_ppe;

% physical velocity
u = real(ifft2(ky .* psi_hat));
v = real(ifft2(-kx .* psi_hat));
dwdx = real(ifft2(kx .* w_hat));
dwdy = real(ifft2(ky .* w_hat));

% convection
% computed in physical space
H = u.*dwdx + v.*dwdy;
% transformed to FFT space
H_hat = fft2(H);
H_hat = cutoff_filter .* H_hat;

% diffusion
D_hat = nu * k2_visc .* w_hat;

dw_hat = H_hat + D_hat;

return
end

