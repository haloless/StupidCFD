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

% ## SpaceTimeDGInit1D

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-01-22

function [ ret ] = SpaceTimeDGInit1D (nvar, use_modal)

% ADER-WENO is needed a priori
ADERWENOGlobals1D;
% to be initialized
SpaceTimeDGGlobals1D;

M = MDegree;
N = NPoint;
Nd = N * N;

eta = GausEta;
wgt = GausWgt;

% TODO NVAR>1
if (nvar~=1)
    error('NVAR=1 only');
end

UseModal = use_modal;

MMat = zeros(Nd,Nd);
KtMat = zeros(Nd,Nd);
KxMat = zeros(Nd,Nd);
F0Mat = zeros(Nd,N);
F1Mat = zeros(Nd,Nd);
tau0 = 0;
tau1 = 1;

K1Mat = zeros(Nd,Nd);

if (use_modal) % modal form
    % error('modal form not available')
    warning('Modal form is only for test');
    
    % mass matrix
    for k = 1:Nd
    for l = 1:Nd
        [ki,kj] = ind2sub([N N], k);
        [li,lj] = ind2sub([N N], l);
        if (ki~=li || kj~=lj); continue; end
        vki = RescaledLegendrePolyVal(ki-1,eta);
        vkj = RescaledLegendrePolyVal(kj-1,eta);
        vli = RescaledLegendrePolyVal(li-1,eta);
        vlj = RescaledLegendrePolyVal(lj-1,eta);
        vi = sum(vki .* vli .* wgt);
        vj = sum(vkj .* vlj .* wgt);
        MMat(k,l) = vi * vj;
    end
    end

    % stiffness matrix
    for k = 1:Nd
    for l = 1:Nd
        [ki,kj] = ind2sub([N N], k);
        [li,lj] = ind2sub([N N], l);
        cki = RescaledLegendrePolyCoef(ki-1);
        ckj = RescaledLegendrePolyCoef(kj-1);
        cli = RescaledLegendrePolyCoef(li-1);
        clj = RescaledLegendrePolyCoef(lj-1);
        vki = polyval(cki,eta);
        vkj = polyval(ckj,eta);
        vli = polyval(cli,eta);
        vlj = polyval(clj,eta);
        dkj = polyval(polyder(ckj),eta);
        dli = polyval(polyder(cli),eta);
        
        KtMat(k,l) = sum(vki.*vli.*wgt) * sum(dkj.*vlj.*wgt);
        KxMat(k,l) = sum(vki.*dli.*wgt) * sum(vkj.*vlj.*wgt);
    end
    end

    % flux matrix TAU=0
    for k = 1:Nd
    for m = 1:N
        [ki,kj] = ind2sub([N N], k);
        if (ki~=m); continue; end
        vki = RescaledLegendrePolyVal(ki-1,eta);
        vkj = RescaledLegendrePolyVal(kj-1,tau0);
        vm = RescaledLegendrePolyVal(m-1,eta);
        F0Mat(k,m) = sum(vki.*vkj.*vm.*wgt);
    end
    end

    % flux matrix TAU=1
    for k = 1:Nd
    for l = 1:Nd
        [ki,kj] = ind2sub([N N], k);
        [li,lj] = ind2sub([N N], l);
        if (ki~=li); continue; end
        vki = RescaledLegendrePolyVal(ki-1,eta);
        vkj = RescaledLegendrePolyVal(kj-1,tau1);
        vli = RescaledLegendrePolyVal(li-1,eta);
        vlj = RescaledLegendrePolyVal(lj-1,tau1);
        F1Mat(k,l) = sum(vki.*vkj.*vli.*vlj.*wgt);
    end
    end
else % nodal form
    % mass matrix
    for q = 1:Nd
    for p = 1:Nd
        [qx,qt] = ind2sub([N N], q);
        [px,pt] = ind2sub([N N], p);
        % if (qx~=px || qt~=pt); continue; end;
        if (0)
            vqx = polyval(LagrPsi(:,qx),eta);
            vqt = polyval(LagrPsi(:,qt),eta);
            vpx = polyval(LagrPsi(:,px),eta);
            vpt = polyval(LagrPsi(:,pt),eta);
            MMat(q,p) = sum(vqx.*vpx.*wgt) * sum(vqt.*vpt.*wgt);
        else
            if (qx==px && qt==pt)
                MMat(q,p) = wgt(px) * wgt(pt);
            end
        end
    end
    end
    % space stiffness matrix
    for q = 1:Nd
    for p = 1:Nd
        [qx,qt] = ind2sub([N N], q);
        [px,pt] = ind2sub([N N], p);
        if (0)
            vqx = polyval(LagrPsi(:,qx),eta);
            vqt = polyval(LagrPsi(:,qt),eta);
            cpx = LagrPsi(:,px);
            dpx = polyval(polyder(cpx),eta);
            vpt = polyval(LagrPsi(:,pt),eta);
            KxMat(q,p) = sum(vqx.*dpx.*wgt) * sum(vqt.*vpt.*wgt);
        else
            if (qt==pt)
                dpx = polyval(polyder(LagrPsi(:,px)),eta(qx));
                KxMat(q,p) = dpx*wgt(qx) * wgt(pt);
            end
        end
    end
    end
    % time stiffness matrix
    for q = 1:Nd
    for p = 1:Nd
        [qx,qt] = ind2sub([N N], q);
        [px,pt] = ind2sub([N N], p);
        if (0)
            vqx = polyval(LagrPsi(:,qx),eta);
            cqt = LagrPsi(:,qt);
            dqt = polyval(polyder(cqt),eta);
            vpx = polyval(LagrPsi(:,px),eta);
            vpt = polyval(LagrPsi(:,pt),eta);
            KtMat(q,p) = sum(vqx.*vpx.*wgt) * sum(dqt.*vpt.*wgt);
        else
            if (qx == px)
                dqt = polyval(polyder(LagrPsi(:,qt)),eta(pt));
                KtMat(q,p) = wgt(px) * dqt*wgt(pt);
            end
        end
    end
    end
    % flux matrix TAU=0
    for q = 1:Nd
    for m = 1:N
        [qx,qt] = ind2sub([N N], q);
        if (0)
            vqx = polyval(LagrPsi(:,qx),eta);
            vqt = polyval(LagrPsi(:,qt),tau0);
            vm = polyval(LagrPsi(:,m),eta);
            F0Mat(q,m) = sum(vqx.*vqt.*vm.*wgt);
        else
            if (qx == m)
                vqt0 = polyval(LagrPsi(:,qt),tau0);
                F0Mat(q,m) = wgt(m) * vqt0;
            end
        end
    end
    end
    % flux matrix TAU=1
    for q = 1:Nd
    for p = 1:Nd
        [qx,qt] = ind2sub([N N], q);
        [px,pt] = ind2sub([N N], p);
        if (0)
            vqx = polyval(LagrPsi(:,qx),eta);
            vqt = polyval(LagrPsi(:,qt),tau1);
            vpx = polyval(LagrPsi(:,px),eta);
            vpt = polyval(LagrPsi(:,pt),tau1);
            F1Mat(q,p) = sum(vqx.*vpx.*vqt.*vpt.*wgt);
        else
            if (qx == px)
                vqt1 = polyval(LagrPsi(:,qt),tau1);
                vpt1 = polyval(LagrPsi(:,pt),tau1);
                F1Mat(q,p) = wgt(px) * vqt1*vpt1;
            end
        end
    end
    end
end

K1Mat = F1Mat - KtMat;

return
end
