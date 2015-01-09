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

% ## TransportTest1D

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-09-06

clc;
clear all;

n = 160;
c = 1.1;
dt = 4e-3;

xlo = -1;
xhi = 1;
xlen = xhi - xlo;
h = xlen / (n-1);
x = linspace(xlo,xhi,n)';

%
diffl = @(x) [0; diff(x)];
diffr = @(x) [diff(x); 0];
sumr = @(x) [x(1:end-1)+x(2:end); x(end)*2];
avgr = @(x) [1/2*(x(1:end-1)+x(2:end)); x(end)];
%
theta_func = @(u) diffl(u) ./ (diffr(u) + eps);
limiter_vonleer = @(r) (abs(r)+r) ./ (1+abs(r));
limiter_superbee = @(r) max(0, max(min(1,2*r), min(r,2)));


%
init_cond_func = @(x) (x>-0.7&x<-0.2)*0.8 + (x>-0.9&x<-0.6)*0.2;
u0 = init_cond_func(x);
% u = u0;
u_upwind = u0;
u_lf = u0;
u_lw = u0;
u_lw_vonleer = u0;
u_lw_superbee = u0;

step = 0;
time = 0;
max_time = 1.0;
while (time<max_time)
    time = time + dt;
    step = step + 1;
    
    % upwind
    f_upwind = c * u_upwind;
    % LF
    f_lf = c*avgr(u_lf) - h/(2*dt)*diffr(u_lf);
    % LW
    f_lw = c*avgr(u_lw) - dt/2*c^2*diffr(u_lw);
    % LW + von Leer
    theta = theta_func(u_lw_vonleer);
    phi = limiter_vonleer(theta);
    f_lw_vonleer = c*u_lw_vonleer + c/2*(1-c*dt/h)*diffr(u_lw_vonleer).*phi;
    % LW + Superbee
    theta = theta_func(u_lw_superbee);
    phi = limiter_superbee(theta);
    f_lw_superbee = c*u_lw_superbee + c/2*(1-c*dt/h)*diffr(u_lw_superbee).*phi;
    
    u_upwind = u_upwind - dt/h*diffl(f_upwind);
    u_lf = u_lf - dt/h*diffl(f_lf);
    u_lw = u_lw - dt/h*diffl(f_lw);
    u_lw_vonleer = u_lw_vonleer - dt/h*diffl(f_lw_vonleer);
    u_lw_superbee = u_lw_superbee - dt/h*diffl(f_lw_superbee);
    
    
    if (1)
        clf;
        plot(x,u0,':', x,init_cond_func(x-c*time),'-', ...
        x,u_upwind,'+-', x,u_lf,'o-', x,u_lw,'*-', ...
        x,u_lw_vonleer,'x-', x,u_lw_superbee,'s-');
        legend('init','trans', ...
        'upwind', 'Lax-Friedrichs','Lax-Wendroff','LW+von Leer','LW+Superbee');
        axis([xlo xhi -0.2 1.2]);
        title(sprintf('step=%d; time=%f',step,time));
        drawnow;
    end
end








