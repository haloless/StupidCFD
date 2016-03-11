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

% ## TestMG1D

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-01-25

% clc
clear all

SemiMGGlobals;
%
smg_verbose = 1;
smg_debug = 1;
smg_nu1 = 1;
smg_nu2 = 1;
smg_nu0 = 1;
%
smg_num_level = -1;

%
n = 8;
% n = 10;
% n = 15;
nx = n;
ny = n;
Lx = 1.0;
Ly = 1.0;
dx = Lx / nx;
dy = Ly / ny;

gravy = -1;


% setup test problem
[ Ac Aw Ae As An rhs ] = TestInitProb(nx,ny,dx,dy);

% zero solution
sol = zeros(nx,ny);


% init solver
nlevel = SemiMGInit(nx,ny);
disp(['SMG nlevel=', int2str(nlevel)]);

% generate coefficients
SemiMGSetup(Ac, Aw, Ae, As, An);


if (smg_debug>0 && nx==8 && ny==8)
%
b1 = rhs;
x1 = sol;

%
x1 = SemiMGSmooth(1, x1, b1);
r1 = SemiMGResid(1, x1, b1);

%
b2 = SemiMGRestrict(2, r1);
x2 = zeros(smg_levels(2).ncell);
x2 = SemiMGSmooth(2, x2, b2);
r2 = SemiMGResid(2, x2, b2);

%
b3 = SemiMGRestrict(3, r2);
x3 = zeros(smg_levels(3).ncell);
x3 = SemiMGSmooth(3, x3, b3);
r3 = SemiMGResid(3, x3, b3);

%
b4 = SemiMGRestrict(4, r3);
x4 = zeros(smg_levels(4).ncell);
x4 = SemiMGSmooth(4, x4, b4);
r4 = SemiMGResid(4, x4, b4);

%
b5 = SemiMGRestrict(5, r4);
x5 = zeros(smg_levels(5).ncell);
x5 = SemiMGSmooth(5, x5, b5);
r5 = SemiMGResid(5, x5, b5);

%
b6 = SemiMGRestrict(6, r5);
x6 = zeros(smg_levels(6).ncell);
x6 = SemiMGSmooth(6, x6, b6);
r6 = SemiMGResid(6, x6, b6);

%
b7 = SemiMGRestrict(7, r6);
x7 = zeros(smg_levels(7).ncell);
x7 = SemiMGSmooth(7, x7, b7);

% 
e6 = SemiMGProlong(6, x7);
z6 = x6 + e6;
x6 = SemiMGSmooth(6, z6, b6);

%
e5 = SemiMGProlong(5, x6);
z5 = x5 + e5;
x5 = SemiMGSmooth(5, z5, b5);

%
e4 = SemiMGProlong(4, x5);
z4 = x4 + e4;
x4 = SemiMGSmooth(4, z4, b4);

%
e3 = SemiMGProlong(3, x4);
z3 = x3 + e3;
x3 = SemiMGSmooth(3, z3, b3);

%
e2 = SemiMGProlong(2, x3);
z2 = x2 + e2;
x2 = SemiMGSmooth(2, z2, b2);

%
e1 = SemiMGProlong(1, x2);
z1 = x1 + e1;
x1 = SemiMGSmooth(1, z1, b1);

else
sol = SemiMGRelax(1, sol, rhs);
end















