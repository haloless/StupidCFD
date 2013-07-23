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

% ## CWENO1D_driver

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-15

function CWENO1D_driver
cfl = 0.4;
L = 2.0;
ncell = 1280;
dx = L / ncell;

ngrow = 4;


% storage
nmax = 2800;
rab = zeros(3,1);
a = zeros(3,3);
b = zeros(3,3);
rf = zeros(3,1);
up = zeros(nmax,3);
wf = zeros(nmax,3,4);
w = zeros(3,1);
s = zeros(3,1);
wab = zeros(nmax,3);
cf = zeros(3,3);
pp1 = zeros(3,1);
pp2 = zeros(3,1);
rc = zeros(3,1);
u = zeros(nmax,1);
rcn = zeros(3,1);
rcp = zeros(3,1);
wcp = zeros(nmax,3);
ff = zeros(nmax,4);
c = zeros(3,3);
v = zeros(nmax,1);
gi = zeros(nmax,4);
wcn = zeros(nmax,3);
ffd = zeros(nmax,1);

% 3-point Gauss, on [0,1]
gauss3_a = [5/18; 8/18; 5/18];
gauss3_x = 0.5 + [-sqrt(15)/10, 0, +sqrt(15)/10];

% RK coefficient
a21 = 0.5;
a32 = 0.5;
a43 = 1.0;
b1 = 1/6;
b2 = 1/3;
b3 = 1/3;
b3 = 1/6;
c2 = 0.5;
c3 = 0.5;
c4 = 1.0;

% coefficient of p_j(x) at [x_j+1/2,x_j+1]
a = [ ...
1/16, -1/4, 11/16; ...
-1/16, 1/2, 1/16; ...
5/16, 1/4, -1/16];
% coefficient of p_j(x) at [x_j,x_j+1/2]
b = [ ...
-1/16, 1/4, 5/16; ...
1/16, 1/2, -1/16; ...
11/16, -1/4, 1/16];

% coefficient of p(x) at point x_j
c = [ ...
-1/24, 1/12, 23/24; ...
-1/24, 13/12, -1/24; ...
23/24, 1/12, -1/24];
% coefficient of p'(x) for f at point x_j
cf = [ ...
1/2, -2, 3/2; ...
-1/2, 0, 1/2; ...
-3/2, 2, -1/2];

end
fprintf('\n')

% ================================================

function [ u0 ] = u0 (x)
u0 = 3/8*x - 1/(4*pi)*sin(2*pi*x) + 1/(32*pi)*sin(4*pi*x);








