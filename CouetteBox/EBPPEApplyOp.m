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

% ## EBPPEApplyOp

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-09-05

function [ y ] = EBPPEApplyOp (nx,ny, Be,Bw,Bn,Bs,Bp, phi)
% Description
% A * phi = y

% check
if (size(phi,1)~=nx+2 || size(phi,2)~=ny+2)
    error('Input wrong size');
end

% apply BC
% FIXME
phi = PressureBC(phi,nx,ny);

I = 2:nx+1;
J = 2:ny+1;

y = zeros(nx+2,ny+2);
y(I,J) = Bp.*phi(I,J) ...
- Be.*phi(I+1,J) - Bw.*phi(I-1,J) ...
- Bn.*phi(I,J+1) - Bs.*phi(I,J-1);
% y(I,J) = Be .* (phi(I,J) - phi(I+1,J)) ...
% + Bw .* (phi(I,J) - phi(I-1,J)) ...
% + Bn .* (phi(I,J) - phi(I,J+1)) ...
% + Bs .* (phi(I,J) - phi(I,J-1));
return
end
