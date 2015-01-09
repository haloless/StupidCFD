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

% ## SWE1D_SystemMatrix

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-02-02

function [ A,D,R,RInv ] = SWE1D_SystemMatrix (ucons)

SWE1D_Globals;

h = ucons(UH);
u = ucons(UHU) / h;
u2 = u^2;
gh = GRAV * h;

A = [0, 1, 0; ...
-u2+gh, u*2, gh; ...
0, 0, 0];

if (nargout >= 2)
    c = sqrt(gh);
    D = diag([u-c, u+c, 0]);
end
if (nargout >= 3)
    R = [1, 1, gh; ...
    u-c, u+c, 0; ...
    0, 0, u2-gh];
    RInv = inv(R);
    
    if(1)
        Ar = A - R*D*RInv;
        if (max(Ar(:)) > 1e-8)
            A
            R*D*RInv
            error('system matrix invalid');
        end
    end
end

return
end
