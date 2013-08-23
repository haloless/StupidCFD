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

% ## ../VorticityBox/TimeIntegrators

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-16

function [ ret ] = TimeIntegrators ()


return
end

% SSPRK(2,2)
function [ un ] = SSPRK_2_2(u,L,dt)
    u1 = u + dt*L(u);
    un = 1/2*u + 1/2*u1 + 1/2*dt*L(u1);
end

% SSPRK(3,3)
function [ un ] = SSPRK_3_3(u,L,dt)
    u1 = u + dt*L(u);
    u2 = 3/4*u + 1/4*u1 + 1/4*dt*L(u1);
    un = 1/3*u + 2/3*u2 + 2/3*dt*L(u2);
return
end

% SSPRK(5,4)
function [ un ] = SSPRK_5_4(u,L,dt)
    u1 = u + 0.391752226571890*dt*L(u);
    u2 = 0.444370493651235*u + 0.555629506348765*u1 + 0.368410593050371*dt*L(u1);
    u3 = 0.620101851488403*u + 0.379898148511597*u2 + 0.251891774271694*dt*L(u2);
    u4 = 0.178079954393132*u + 0.821920045606868*u3 + 0.544974750228521*dt*L(u3);
    un = 0.517231671970585*u2 ...
        + 0.096059710526147*u3 + 0.063692468666290*dt*L(u3) ...
        + 0.386708617503269*u4 + 0.226007483236906*dt*L(u4);
return
end




