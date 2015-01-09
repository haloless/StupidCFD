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

% ## EBPPEGenerateCoef

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-09-05

function [ Be,Bw,Bn,Bs,Bp ] = EBPPEGenerateCoef (nx,ny,dx,dy, adiag,bx,by)
% Description

% check
if (size(adiag,1)~=nx || size(adiag,2)~=ny)
    error('Alpha wrong size');
end
if (size(bx,1)~=nx+1 || size(bx,2)~=ny)
    error('Beta_x wrong size');
end
if (size(by,1)~=nx || size(by,2)~=ny+1)
    error('Beta_y wrong size');
end

% generate coefficients
I = 1:nx;
J = 1:ny;
Be = 1/dx^2 * bx(I+1,J);
Bw = 1/dx^2 * bx(I,J);
Bn = 1/dy^2 * by(I,J+1);
Bs = 1/dy^2 * by(I,J);
Bp = adiag + (Be + Bw + Bn + Bs);

return
end
