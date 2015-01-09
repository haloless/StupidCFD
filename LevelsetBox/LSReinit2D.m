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

% ## LSReinit2D

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-04-27

function [ dd ] = LSReinit2D (dd,lo,hi,ng,dh)

LSGlobals2D;

% assume input LS is properly filled at the ghost cells
% now D4 is a copy of DD
dd = LSFillGhost2D(dd,lo,hi,ng);
d4 = dd;

% TODO get area

nband = dens_spread + 1;
error_width = nband * dh;
tau = 0.0; % artificial time
tau_max = error_width;
dtau = dh / dhdiv;
step_limit = floor(nband * dhdiv);
max_steps = -1;
% max_steps = 1;
% max_steps = 16;

error_small = 0;
steps = 0;
while (~error_small)
    tau = tau + dtau;
    steps = steps + 1;
    disp(['Re-initialization:step=', int2str(steps), ';tau/dh=', num2str(tau/dh)]);
    
    % fill ghost cells for D4
    d4 = LSFillGhost2D(d4,lo,hi,ng);
    d5 = d4; % D5 holds the initial value for D4 before each RK2 stepping
    
    for rk_count = 0:1 % RK sub-step
        %
        d4 = LSFillGhost2D(d4,lo,hi,ng);
        
        % compute |d|
        [yy,max_radius,d_error,d_count] = LSNormGrad2D(dd,d4,lo,hi,ng,dh,error_width);
        
        % TODO store D6
        
        for j = lo(2):hi(2)
        for i = lo(1):hi(1)
            main_level = dd(i,j);
            curr_level = d4(i,j);
            curr_radius = yy(i,j);
            bias = LSBiasSign(main_level, dh);
            
            dchange = bias * (1.0-curr_radius);
            
            if (rk_count == 0) % 1st RK2 stage
                new_level = curr_level + dchange*dtau;
                % d5(i,j) = curr_level;
            else % 2nd RK2 stage
                new_level = 0.5 * (d5(i,j) + curr_level+dchange*dtau);
            end
            
            if (new_level<0.0 && main_level>=0.0)
                new_level = 0.0;
            elseif (new_level>0.0 && main_level<=0.0)
                new_level = 0.0;
            end
            
            d4(i,j) = new_level;
        end
        end
    end % RK_COUNT
    
    % TODO conserve mass
    
    if (max_steps>0 && steps>=max_steps)
        error_small = 1;
    elseif (tau >= tau_max)
        error_small = 1;
    elseif (steps >= step_limit)
        error_small = 1;
    else
        error_small = 0;
    end
end

% TODO check area

dd = d4;

return
end



