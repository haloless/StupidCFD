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

% ## Burgers1DMain

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-01



% Godunov's method for inviscid Burgers' equation

cfl = 0.9;

xa = 0; xb = 1;
xlen = xb-xa;
ncell = 100;


% initilize
nx = ncell + 2;
dx = xlen / ncell;
flux = zeros(ncell+1,1);
U = zeros(nx,1);
X = linspace(xa-dx/2,xb+dx/2,nx)';

prob = 'SmoothGaussian';
% prob = 'SquaredWave'
bc_type = 0; % periodic

switch prob
case {'SmoothGaussian'}
    U0 = exp(-8 * (2*(X/xlen-0.5)).^2);
case {'SquaredWave'}
    U0 = -1.0 * (X<=0.1*xlen) + ...
        1.0 * (0.1*xlen<X & X<=0.5*xlen) + ...
        0.0 * (0.5*xlen<X & X<=0.9*xlen) - ...
        1.0 * (0.9*xlen<X);
otherwise
    error('Unknown prob: %s',prob);
end

U = U0;

if (1)
    figure; plot(X,U,'x-'); title('initial cond.');
end



% computation begins
time = 0;
step = 0;
max_time = 0.5;
max_step = 100000;
while (time<max_time && step<max_step)
    % apply BC
    U = BoundaryCondition1D(U,ncell,1,bc_type);
    
    % estimate dt
    umax = max(abs(U));
    dt = cfl * dx / umax;
    if (dt+time > max_time)
        dt = max_time - time;
    end
    
    time = time + dt;
    step = step + 1;
    
    % numerical flux at cell face
    for i = 1:ncell+1
        ul = U(i);
        ur = U(i+1);
        ustar = BurgersExactRiemann1D(ul,ur);
        flux(i) = 0.5 * ustar^2;
    end
    
    % update
    U(2:ncell+1) = U(2:ncell+1) - dt/dx * diff(flux);
    
    
    if (mod(step,1) == 0)
        disp(['step=',int2str(step), ';time=',num2str(time)]);
        
        plot(X(2:ncell+1),U(2:ncell+1)); title(num2str(time));
        axis([xa xb min(U0) max(U0)]);
        drawnow;
        pause(0.05);
    end
end






