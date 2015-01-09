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

% ## WENO3_plus

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-08

function [ ue ] = WENO3_upwind (u)

% central stencil
i = 3;

% u should be a column vector
if (isrow(u))
    u = u';
end

C = [...
2/6, -7/6, 11/6; ...
-1/6, 5/6, 2/6; ...
2/6, 5/6, -1/6];


usn = [u(i-2:i), u(i-1:i+1), u(i:i+2)];

% polynomials
psn = zeros(3,1);
for s = 1:3
    psn(s) = C(s,:)*usn(:,s);
end

% smoothness indicator
Ba = [1,-2,1; 1,-2,1; 1,-2,1];
Bb = [1,-4,3; 1,0,-1; 3,-4,1];
Bsn = zeros(3,1);
for s = 1:3
    Bsn(s) = 13/12*(Ba(s,:)*usn(:,s))^2 + 1/4*(Bb(s,:)*usn(:,s))^2;
end

dsn = [1/10; 6/10; 3/10];
eps = 1e-6;
% eps = 1e-10;
% alpha weights
alphasn = dsn ./ (eps + Bsn).^2;

% stencils weights
wsn = alphasn ./ sum(alphasn);

% approximation of u_minus_{1/2}
ue = wsn' * psn;

return
end
