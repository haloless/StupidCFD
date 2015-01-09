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

% ## TwoLayerSWE1D_LuOp

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-02-04

function [ Lu,dtau ] = TwoLayerSWE1D_LuOp (ucons, lo,hi,nx,ng,hx)

SWE1D_Globals;

% path integral
gpnt = [1/2-sqrt(15)/10, 1/2, 1/2+sqrt(15)/10]';
gwgt = [5/18, 8/18, 5/18]';
gnum = 3;

% apply BC
ucons = SWE1D_FillBC(ucons, lo,hi,ng);

% reconstruction
urcl = ucons;
urch = ucons;
urcp = zeros(NRECONS,NUCONS,nx);
% reconstruction correction
Ad = zeros(NUCONS,nx);

if (NRECONS == 1)
    urcp(1,:,:) = ucons;
elseif (NRECONS == 2)
    urcd = zeros(NRECONS-1,NUCONS,nx);
    uslp = zeros(NUCONS,nx);
    for i = lo-1:hi+1
        dupw = ucons(:,i) - ucons(:,i-1);
        dloc = ucons(:,i+1) - ucons(:,i);
        % limited slope
        uslp(:,i) = slope_minmod(dupw,dloc);
    end
    % extrapolate to cell face
    urcl = ucons - 0.5*uslp;
    urch = ucons + 0.5*uslp;
    
    for n = 1:gnum
        uply = urcl + gpnt(n)*uslp;
        for i = lo:hi
            An = TwoLayerSWE1D_SystemMatrix(uply(:,i));
            Ad(:,i) = Ad(:,i) + gwgt(n)*An*uslp(:,i);
        end
    end

    % % save reconstructed polynomial and derivative
    % urcp(1,:,:) = uslp;
    % urcp(2,:,:) = urcl;
    % urcd(1,:,:) = uslp;
end

% if (NRECONS > 1)
    % for i = lo:hi
        % upval = zeros(NUCONS,gnum);
        % udval = zeros(NUCONS,gnum);
        % for comp = 1:NUCONS
            % uply = reshape(urcp(:,comp,i),[],1);
            % uder = reshape(urcd(:,comp,i),[],1);
            % upval(comp,:) = polyval(uply,gpnt);
            % udval(comp,:) = polyval(uder,gpnt);
        % end
        
        % for n = 1:gnum
            % An = TwoLayerSWE1D_SystemMatrix(upval(:,n));
            % Ad(:,i) = Ad(:,i) + gwgt(n)*An*udval(:,n);
        % end
    % end
% end


% jump term
Dp = zeros(NUCONS,nx+1);
Dm = zeros(NUCONS,nx+1);

for i = lo:hi+1 % loop cell face
    ul = urch(:,i-1);
    ur = urcl(:,i);
    
    Ap = zeros(NUCONS,NUCONS);
    Am = zeros(NUCONS,NUCONS);
    for n = 1:gnum % Osher-type
        s = gpnt(n); w = gwgt(n);
        us = ul + s*(ur-ul);
        
        [As,Ds,Rs,Ls] = TwoLayerSWE1D_SystemMatrix(us);
        Aa = Rs * abs(Ds) * Ls;
        
        Ap = Ap + w/2 * (As + Aa);
        Am = Am + w/2 * (As - Aa);
    end
    
    Dp(:,i) = Ap * (ur-ul);
    Dm(:,i) = Am * (ur-ul);
end


I = lo:hi;
Lu = zeros(size(ucons));
Lu(:,I) = -1/hx * (Dm(:,I+1) + Dp(:,I) + Ad(:,I));

% estimate time step
dtau = TwoLayerSWE1D_EstimTimeStep(ucons,lo,hi,hx);

return
end

function [ dtau ] = TwoLayerSWE1D_EstimTimeStep(ucons, lo,hi,hx)
SWE1D_Globals;
%
c1s = sqrt(GRAV .* ucons(UH1,:));
u1s = ucons(UHU1,:) ./ ucons(UH1,:);
a1s = abs(u1s) + c1s;
%
c2s = sqrt(GRAV .* ucons(UH2,:));
u2s = ucons(UHU2,:) ./ ucons(UH2,:);
a2s = abs(u2s) + c2s;
%
dtau = hx / max(max(a1s),max(a2s));

return
end % TwoLayerSWE1D_EstimTimeStep

function [ A,D,R,L ] = TwoLayerSWE1D_SystemMatrix(s)
SWE1D_Globals;
h1 = s(UH1); u1 = s(UHU1) / h1;
h2 = s(UH2); u2 = s(UHU2) / h2;
gh1 = GRAV * h1;
gh2 = GRAV * h2;

A = [0, 1, 0, 0, 0; ...
-u1^2+gh1, 2*u1, gh1, 0, gh1; ...
0, 0, 0, 1, 0; ...
DENS*gh2, 0, -u2^2+gh2, 2*u2, gh2; ...
0, 0, 0, 0, 0];

if (nargout > 1)
    [R,D] = eig(A);
    L = inv(R);
    if (1)
        A = real(A); D = real(D); R = real(R); L = real(L);
    end
end
return
end % TwoLayerSWE1D_SystemMatrix

function [s] = slope_minmod(sl,sr)
% s = zeros(size(sl));
s = (sl.*sr>0) .* sign(sl+sr) .* min(abs(sl),abs(sr));
return
end % slope_minmod


