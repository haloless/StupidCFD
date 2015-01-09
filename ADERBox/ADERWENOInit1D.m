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

% ## ADERWENOInit1D

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-01-19

function [ ret ] = ADERWENOInit1D (ader_order)
% Description

ADERWENOGlobals1D;

% order of ADER
N = ader_order;
M = ader_order - 1;
NPoint = N;
MDegree = M;
NGrow = M;

% number of stencils and stencil offsets
lambda0 = 1.0e5;
lambda1 = 1.0;
if (mod(N,2) == 1)
    Ns = 3; % odd order
    soff = [-M, -M/2, 0];
    slambda = [lambda1, lambda0, lambda1];
else
    Ns = 4; % even order
    soff = [-M, -N/2, -N/2+1, 0];
    slambda = [lambda1, lambda0, lambda0, lambda1];
end
NStencil = Ns;
StencilOffset = zeros(N,Ns);
StencilIndice = zeros(N,Ns);
for s = 1:Ns
    StencilOffset(:,s) = (0:M) + soff(s);
    StencilIndice(:,s) = N + StencilOffset(:,s);
end
StencilLambda = slambda';
StencilEpsilon = 1.0e-14;
StencilRaise = 8;
% StencilRaise = 12;

%
[GausEta,GausWgt] = GaussQuadCoefs1D(N, 0, 1);
%
LagrPsi = LagInterpPoly(GausEta);
LagrPsi0 = zeros(N,1);
LagrPsi1 = zeros(N,1);
for i = 1:N
    LagrPsi0(i) = polyval(LagrPsi(:,i),0.0);
    LagrPsi1(i) = polyval(LagrPsi(:,i),1.0);
end
% derivative of Lagr. Interp. Poly.
LagrPsiDer1 = zeros(N-1,N);
for i = 1:N
    LagrPsiDer1(:,i) = polyder(LagrPsi(:,i));
end

%
LagrNodalMat = zeros(N,N,Ns);
LagrNodalMatInv = zeros(N,N,Ns);
for s = 1:Ns
    smat = zeros(N,N);
    for e = 1:N % cells in current stencil
        xe = GausEta + StencilOffset(e,s);
        ce = zeros(N,N);
        for i = 1:N % points in current cell
            for j = 1:N
                ce(i,j) = polyval(LagrPsi(:,j), xe(i));
            end
        end
        smat(e,:) = GausWgt' * ce;
    end
    LagrNodalMat(:,:,s) = smat;
    LagrNodalMatInv(:,:,s) = inv(smat);
end

%
SigmaMat = zeros(N,N);
derlagp = LagrPsi;
for a = 1:M
    derlagp1 = zeros(N-a,N);
    for i = 1:N
        derlagp1(:,i) = polyder(derlagp(:,i));
    end
    derlagp = derlagp1;
    for i = 1:N
    for j = 1:N
        dpi = polyval(derlagp(:,i),GausEta);
        dpj = polyval(derlagp(:,j),GausEta);
        SigmaMat(i,j) = SigmaMat(i,j) + sum(dpi.*dpj.*GausWgt);
    end
    end
end

return
end
