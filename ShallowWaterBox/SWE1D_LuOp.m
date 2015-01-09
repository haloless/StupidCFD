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

% ## SWE1D_LuOp

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-02-02

function [ Lu,dtau ] = SWE1D_LuOp (ucons, lo,hi,nx,ng,hx)

SWE1D_Globals;

% apply BC
ucons = SWE1D_FillBC(ucons, lo,hi,ng);
% get primitive variables
uprim = SWE1D_ConsToPrim(ucons);

% reconstruction
urcl = ucons;
urch = ucons;
urcp = zeros(NRECONS,NUCONS,nx);
if (NRECONS == 1)
    urcp(1,:,:) = ucons;
elseif (NRECONS == 2) % 2nd order TVD
    urcd = zeros(NRECONS-1,NUCONS,nx);
    for i = lo-1:hi+1 % loop cell with 1 ghost
        dupw = ucons(:,i) - ucons(:,i-1);
        dloc = ucons(:,i+1) - ucons(:,i);
        delta = 0.5*dupw + 0.5*dloc;
        % limited slope
        for comp = 1:NUCONS
            delta(comp) = slope_minmod(dupw(comp),dloc(comp));
            % delta(comp) = slope_mc(dupw(comp),dloc(comp));
            % delta(comp) = slope_superbee(dupw(comp),dloc(comp));
        end
        % extrapolate to cell face
        urcl(:,i) = ucons(:,i) - 0.5*delta;
        urch(:,i) = ucons(:,i) + 0.5*delta;
        % save reconstructed polynomial
        urcp(1,:,i) = delta;
        urcp(2,:,i) = urcl(:,i);
        % derivative
        urcd(1,:,i) = delta;
    end
elseif (NRECONS >= 3) % ADER-WENO
    urcd = zeros(NRECONS-1,NUCONS,nx);
    M = NRECONS - 1;
    for i = lo-1:hi+1
        uloc = ucons(:,i-M:i+M);
        if (0) % component-wise 
            for comp = 1:NUCONS
                % [qw,pw] = ADERWENOReconstruct1D(ucons(comp,i-M:i+M)',hx);
                [qw,pw] = ADERWENOReconstruct1D(uloc(comp,:)',hx);
                %
                urcl(comp,i) = polyval(pw,0);
                urch(comp,i) = polyval(pw,1);
                %
                urcp(:,comp,i) = pw;
                urcd(:,comp,i) = polyder(pw);
            end
        else % characteristic-wise
            [As,Ds,Rs,Ls] = SWE1D_SystemMatrix(ucons(:,i));
            uloc = Ls * uloc;
            uqw = zeros(NUCONS,NRECONS);
            for comp = 1:NUCONS
                [qw,pw] = ADERWENOReconstruct1D(uloc(comp,:)',hx);
                
                uqw(comp,:) = qw;
                % %
                % urcl(comp,i) = polyval(pw,0);
                % urch(comp,i) = polyval(pw,1);
                % %
                % urcp(:,comp,i) = pw;
                % urcd(:,comp,i) = polyder(pw);
            end
            % transform
            uqw = Rs * uqw;
            for comp = 1:NUCONS
                upw = ADERWENOLagrPoly1D(uqw(comp,:)', 0);
                udw = ADERWENOLagrPoly1D(uqw(comp,:)', 1);
                %
                urcl(comp,i) = polyval(upw,0);
                urch(comp,i) = polyval(upw,1);
                %
                urcp(:,comp,i) = upw;
                urcd(:,comp,i) = udw;
            end
        end
    end
else % invalid reconstruction
    error('Invalid reconstruction');
end

% path integral
if (1)
    gpnt = [1/2-sqrt(15)/10, 1/2, 1/2+sqrt(15)/10]';
    gwgt = [5/18, 8/18, 5/18]';
    gnum = 3;
else
    gpnt = [-1/21*sqrt(245+14*sqrt(70)), -1/21*sqrt(245-14*sqrt(70)), 0, ...
    1/21*sqrt(245-14*sqrt(70)), 1/21*sqrt(245+14*sqrt(70))]';
    gwgt = [1/900*(322-13*sqrt(70)), 1/900*(322+13*sqrt(70)), 128/225, ...
    1/900*(322+13*sqrt(70)), 1/900*(322-13*sqrt(70))];
    gpnt = 0.5*gpnt + 0.5;
    gwgt = 0.5*gwgt;
    gnum = 5;
end


% jump term
Dp = zeros(NUCONS,nx+1);
Dm = zeros(NUCONS,nx+1);

for i = lo:hi+1 % loop cell face
    ul = urch(:,i-1);
    ur = urcl(:,i);
    
    Ap = zeros(NUCONS,NUCONS);
    Am = zeros(NUCONS,NUCONS);
    for n = 1:gnum % Osher-type Riemann solver
        s = gpnt(n);
        w = gwgt(n);
        us = ul + s*(ur-ul);
        [As,Ds,Rs,RsInv] = SWE1D_SystemMatrix(us);
        Aa = Rs * abs(Ds) * RsInv;
        
        Ap = Ap + w/2 * (As + Aa);
        Am = Am + w/2 * (As - Aa);
    end
    if (0) % switch to Roe-type scheme
        As = Ap + Am;
        [Rs,Ds] = eig(As);
        Aa = Rs * abs(Ds) * inv(Rs);
        Ap = 1/2 * (As + Aa);
        Am = 1/2 * (As - Aa);
    end
    
    Dp(:,i) = Ap * (ur-ul);
    Dm(:,i) = Am * (ur-ul);
end

% reconstruction correction
Ad = zeros(NUCONS,nx);
if(NRECONS > 1)
    for i = lo:hi
        upval = zeros(NUCONS,gnum);
        udval = zeros(NUCONS,gnum);
        for comp = 1:NUCONS
            uply = reshape(urcp(:,comp,i),[],1);
            uder = reshape(urcd(:,comp,i),[],1);
            
            upval(comp,:) = polyval(uply,gpnt);
            udval(comp,:) = polyval(uder,gpnt);
        end
        
        for n = 1:gnum
            An = SWE1D_SystemMatrix(upval(:,n));
            Ad(:,i) = Ad(:,i) + gwgt(n) * An * udval(:,n);
        end
    end
end

I = lo:hi;
Lu = zeros(size(ucons));
Lu(:,I) = -1/hx * (Dm(:,I+1) + Dp(:,I) + Ad(:,I));
% for i = lo:hi
    % Lu(:,i) = -1/hx * (Dm(:,i+1) + Dp(:,i));
% end

asig = abs(uprim(QVX,:)) + uprim(QC,:);
dtau = hx / (max(asig));

return
end


function [R,D,L] = SWE1D_RoeMatrix(sl,sr)
SWE1D_Globals;

hl = sl(UH); 
hr = sr(UH);
ul = sl(UHU) / hl;
ur = sr(UHU) / hr;

wl = sqrt(hl);
wr = sqrt(hr);

hbar = 0.5 * (hl+hr);
ubar = (wl*ul+wr*ur) / (wl+wr);
cbar = sqrt(GRAV * hbar);

D = diag([ubar-cbar, ubar+cbar, 0]);
R = [1, 1, 0; ...
ubar-cbar, ubar+cbar, 0; ...
0, 0, 1];
L = inv(R);

return
end % SWE1D_RoeMatrix

function [s] = slope_minmod(sl,sr)
if (sl*sr > 0)
    s = sign(sl+sr) * min(abs(sl),abs(sr));
else
    s = 0;
end
return
end % slope_minmod
function [s] = slope_mc(sl,sr)
if (sl*sr > 0)
    s = sign(sl+sr) * min([2*abs(sl), 2*abs(sr), 0.5*abs(sl+sr)]);
else
    s = 0;
end
return
end % slope_mc
function [s] = slope_superbee(sl,sr)
if (sl*sr > 0)
    s = sign(sl+sr) * max(min(2*abs(sl),abs(sr)),min(abs(sl),2*abs(sr)));
else
    s = 0;
end
return
end % slope_superbee

