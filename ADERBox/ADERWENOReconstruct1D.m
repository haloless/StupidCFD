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

% ## ADERWENOReconstruct1D

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-01-19

function [ qweno,pweno, qs,ps,sigmas,omegas,ws ] = ADERWENOReconstruct1D (fs, hx)
% Description

ADERWENOGlobals1D;

if (isrow(fs)); fs = fs'; end

% reconstruct polynomials in each stencil
qs = zeros(NPoint,NStencil); % point value
ps = zeros(NPoint,NStencil); % polynomial
sigmas = zeros(NStencil,1); % oscillation indicator
for s = 1:NStencil
    qval = LagrNodalMatInv(:,:,s) * fs(StencilIndice(:,s));
    qs(:,s) = qval;
    ps(:,s) = SumPoly(LagrPsi, qval);
    sigmas(s) = qval' * SigmaMat * qval;
end

% compute stencil weights
omegas = StencilLambda ./ ((sigmas+StencilEpsilon).^StencilRaise);
ws = omegas ./ sum(omegas);

% WENO reconstruction
qweno = qs * ws;
pweno = SumPoly(ps,ws);


return
end
