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

% ## ShallowWater1D_Globals

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-02-02

% gravity
global GRAV
global HEIGHT_SMALL

% BC
global bctype bcfill
global BC_PER BC_NEU BC_DIR

% conservative
global UH UHU UBOT
global NUCONS
% primitive
global QH QVX QBOT QGH QC 
global NQPRIM

% reconstruction order
global NRECONS

% two-layer SWE
global UH1 UHU1 UH2 UHU2 UBOT
global QH1 QVX1 QH2 QVX2 QBOT QGH1 QC1 QGH2 QC2
% density ratio
global DENS
