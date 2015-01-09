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

% ## Euler1D_RoeMatrix

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-01-03

function [ R,D,L ] = Euler1D_RoeMatrix (uconsl,uconsr, upriml,uprimr)

Euler1D_globals;

rhol = uconsl(URHO);
rhor = uconsr(URHO);

vl = upriml(QVX);
vr = uprimr(QVX);

% total specific enthalpy
hl = (uconsl(UETOT) + upriml(QPRES)) / rhol;
hr = (uconsr(UETOT) + uprimr(QPRES)) / rhor;

wl = sqrt(rhol);
wr = sqrt(rhor);

% mean state
vbar = (wl*vl + wr*vr) / (wl+wr);
hbar = (wl*hl + wr*hr) / (wl+wr);
ekbar = 0.5 * vbar^2;
eibar = hbar - ekbar;
abar = sqrt((GAMMA-1) * (hbar-0.5*vbar^2));

D = [vbar-abar, vbar, vbar+abar];
R = [ ...
1.0,            1.0,            1.0;
vbar-abar,      vbar,           vbar+abar;
hbar-vbar*abar, 0.5*vbar^2,     hbar+vbar*abar];
% L = inv(R);
ih = 1 / (hbar - 0.5*vbar^2);
L = [ ...
vbar/(2*abar)+0.25*vbar^2*ih, -1/(2*abar)-0.5*vbar*ih, 0.5*ih;
1-0.5*ih*vbar^2, ih*vbar, -ih;
-vbar/(2*abar)+0.25*vbar^2*ih, 1/(2*abar)-0.5*vbar*ih, 0.5*ih];

% check
if (0)
    if (norm(L*R-eye(3)) > 1e-12)
        error('inverse of right eigenmatrix invalid')
    end
end

return
end

