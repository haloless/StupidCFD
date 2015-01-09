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

% ## LBMMain_LidDrivenCavity

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-13

function [velx,vely] = IBLBM_CavityMain
% clc;
% clear all;

% D2Q9
qwgt = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36];
qex  = [0, 1, 0, -1, 0, 1, -1, -1, 1];
qey  = [0, 0, 1, 0, -1, 1, 1, -1, -1];
qord = [1, 2, 3, 4, 5,  6, 7, 8,  9];
qopp = [1, 4, 5, 2, 3,  8, 9, 6,  7];


refine = 1;
Lx = 128 * refine;
Ly = 128 * refine;
nx = Lx;
ny = Ly;

% number of ghost
ng = 4;
Nx = nx + ng*2;
Ny = ny + ng*2;

hx = Lx / nx;
hy = Ly / ny;
hh = min([hx hy]);
qc = 1;
dt = hh / qc;
cs = 1/sqrt(3) * qc;
cs2 = 1/3 * qc^2;

% determine parameters
Re = 400;
% relaxation parameter, must < 2
% omega = 1.7;
% U0 = 0.05;
U0 = 0.1;

if (exist('omega','var'))
    tau = 1 / omega;
    nu = (tau-0.5) * cs2 * dt;
    U0 = Re*nu/Lx;
elseif (exist('U0','var'))
    nu = U0*Lx / Re;
    tau = 1/cs2*nu/dt + 0.5;
    omega = 1/tau;
else
    nu = 1e-3;
    U0 = Re * nu / Lx;
    tau = 1/cs2*nu/dt + 0.5;
    omega = 1/tau;
end

% compute some parameters
Re_h = U0 * hh / nu;
Cfl = dt*U0 / hh;
Dfl = dt*nu / hh^2;
Mach = U0 / cs;

% setup IB
validI = ng+1:Nx-ng;
validJ = ng+1:Ny-ng;
ibMask = ones(Nx,Ny);
ibMask(validI,validJ) = 0;
ibIdx = find(ibMask);

ibFrac = zeros(Nx,Ny); ibFrac(ibIdx) = 1;
ibU = zeros(Nx,Ny); 
% ibU(ng+1:Nx-ng,Ny-ng+1:Ny) = U0;
ibU(:,Ny-ng+1:Ny) = U0;
ibV = zeros(Nx,Ny);
%
ibFrac = reshape(ibFrac, Nx*Ny,1);
ibFracC = 1 - ibFrac;
ibU = reshape(ibU, Nx*Ny,1);
ibV = reshape(ibV, Nx*Ny,1);


% initial condition
rho0 = ones(Nx*Ny,1);
pres0 = cs2 * rho0;

rho = rho0;
pres = pres0;
u = ibFracC.*zeros(Nx*Ny,1) + ibFrac.*ibU; 
v = ibFracC.*zeros(Nx*Ny,1) + ibFrac.*ibV; 

% pressure-based distrib.
ps = zeros(Nx*Ny,9);
% initialize with equilibrium distribution
ps = LBE_pEq(pres,pres0,u,v, qwgt,qex,qey,qc);

max_step = 1000000;
max_time = max_step*dt;
time = 0;
step = 0;
while (step<max_step && time<max_time)
    time = time + dt;
    step = step + 1;
    
    uold = u;
    vold = v;
    
    % predictor
    % collision
    pEq = LBE_pEq(pres,pres0,u,v, qwgt,qex,qey,qc);
    % pStar = ps - omega*(ps - pEq);
    pStar = (1-omega)*ps + omega*pEq;
    % streaming
    pStar = LBE_pStream(pStar, Nx,Ny,qex,qey);
    
    % corrector
    [pres,u,v] = LBEMacroVars(pStar, pres0, qex,qey,qc);
    u = ibFracC.*u + ibFrac.*ibU;
    v = ibFracC.*v + ibFrac.*ibV;
    
    pEq = LBE_pEq(pres,pres0,u,v, qwgt,qex,qey,qc);
    for i = 1:9
        ps(:,i) = ibFracC.*pStar(:,i) + ibFrac.*pEq(:,i);
    end
    
    pres = sum(ps,2);
    % [pres,u,v] = LBEMacroVars(ps, pres0, qex,qey,qc);
    
    % % MICRO BC
    % % fs = CavityMicroBC(fs,idx_xlo,idx_xhi,idx_ylo,idx_yhi,Ubc,Vbc,qc);
    % ps = CavityMicroBC(ps,idx_xlo,idx_xhi,idx_ylo,idx_yhi,Ubc,Vbc,qc);
    
    % % MACRO
    % % [rho,us,vs] = LBEMacroVars(fs, qex,qey,qc);
    % [pres,us,vs] = LBEMacroVars(ps, pres0, qex,qey,qc);
    
    % % MACRO BC
    % % [rho,us,vs] = CavityMacroBC(rho,us,vs, idx_xlo,idx_xhi,idx_ylo,idx_yhi,Ubc,Vbc);
    % [pres,us,vs] = CavityMacroBC(pres,us,vs, idx_xlo,idx_xhi,idx_ylo,idx_yhi,Ubc,Vbc);
    
    % check convergence
    eps_abs = 1e-7 * U0;
    res_abs = norm([u-uold; v-vold],Inf);
    is_conv = res_abs<eps_abs;
    
    if (mod(step,100)==0 || is_conv)
        prompt = ['step=',int2str(step), ';time=',num2str(time), ...
        ';omega=',num2str(omega), ';tau=',num2str(tau), ';Mach=',num2str(Mach), ...
        ';Re=',num2str(Re), ';U0=',num2str(U0), ';nu=',num2str(nu), ...
        ';|res/tol|=',num2str(res_abs/eps_abs), ...
        ];
        disp(prompt);
        
        velx = reshape(u,Nx,Ny);
        vely = reshape(v,Nx,Ny);
        
        velmag = sqrt(velx.^2 + vely.^2);
        contourf(velmag'./U0, 32);
        colorbar; 
        shading flat;
        axis equal;
        title(prompt);
                
        drawnow;
        
        if (is_conv); break; end
    end
end % end of main loop

if (1)
    GhiaU = [0	0
    -0.0825688073	0.0609756098
    -0.1467889908	0.1077235772
    -0.244648318	0.1768292683
    -0.3272171254	0.2845528455
    -0.1743119266	0.4552845528
    0.0183486239	0.6158536585
    0.1590214067	0.7317073171
    0.2874617737	0.8536585366
    0.5565749235	0.9471544715
    0.6177370031	0.9552845528
    0.6819571865	0.9634146341
    0.755351682	0.9715447154
    0.996941896	0.9959349593
    ];
    GhiaV = [0.0081466395	-0.0020366599
    0.0712830957	0.1812627291
    0.1018329939	0.2281059063
    0.1629327902	0.2790224033
    0.2342158859	0.3014256619
    0.5050916497	0.0488798371
    0.8044806517	-0.3890020367
    0.8594704684	-0.450101833
    0.9083503055	-0.3380855397
    0.9429735234	-0.2301425662
    0.9592668024	-0.1568228106
    0.967413442	-0.1242362525
    1	0
    ];
    
    
    velx = velx(validI,validJ);
    vely = vely(validI,validJ);
    
    uprob = 1/U0 * 1/2*(velx(round(nx/2),:) + velx(round((nx+1)/2),:));
    vprob = 1/U0 * 1/2*(vely(:,round(ny/2)) + vely(:,round((ny+1)/2)));
    
    figure;
    plot(uprob,((1:ny)-0.5)/ny,'-', GhiaU(:,1),GhiaU(:,2),'o'); 
    legend('LBM','Ghia'); title('U'); 
    axis equal; axis([-0.5 1 0 1]);
    
    figure;
    plot(((1:nx)-0.5)/nx,vprob,'-', GhiaV(:,1),GhiaV(:,2),'o'); 
    legend('LBM','Ghia'); title('V'); 
    axis equal; axis([0 1 -0.5 0.5])
end

return
end % end of main function


function [fs] = CavityMicroBC(fs,idx_xlo,idx_xhi,idx_ylo,idx_yhi,Ubc,Vbc,qc)
% MICRO BC
% x-low wall, bounce back
fs(idx_xlo,2) = fs(idx_xlo,4);
fs(idx_xlo,6) = fs(idx_xlo,8);
fs(idx_xlo,9) = fs(idx_xlo,7);
% x-high wall, bounce back
fs(idx_xhi,4) = fs(idx_xhi,2);
fs(idx_xhi,8) = fs(idx_xhi,6);
fs(idx_xhi,7) = fs(idx_xhi,9);
% y-low wall, bounce back
fs(idx_ylo,3) = fs(idx_ylo,5);
fs(idx_ylo,6) = fs(idx_ylo,8);
fs(idx_ylo,7) = fs(idx_ylo,9);
% y-high wall, Zou-He
uLid = Ubc(idx_yhi) ./ qc;
vLid = Vbc(idx_yhi) ./ qc;
rhoLid = 1 ./ (1+vLid) .* ...
    (sum(fs(idx_yhi,[1,2,4]),2) + 2*sum(fs(idx_yhi,[3,6,7]),2));
% fs(idx_yhi,5) = fs(idx_yhi,3) - 2/3/qc*rhoLid.*vLid;
% fs(idx_yhi,9) = fs(idx_yhi,7) + 1/2*(fs(idx_yhi,4)-fs(idx_yhi,2)) ...
    % + 1/2/qc*rhoLid.*uLid - 1/6/qc*rhoLid.*vLid;
% fs(idx_yhi,8) = fs(idx_yhi,6) + 1/2*(fs(idx_yhi,2)-fs(idx_yhi,4)) ...
    % - 1/2/qc*rhoLid.*uLid - 1/6/qc*rhoLid.*vLid;
fs(idx_yhi,5) = fs(idx_yhi,3);
fs(idx_yhi,9) = fs(idx_yhi,7) + 1/2*(fs(idx_yhi,4)-fs(idx_yhi,2)) ...
    + 1/2*rhoLid.*uLid;
fs(idx_yhi,8) = fs(idx_yhi,6) + 1/2*(fs(idx_yhi,2)-fs(idx_yhi,4)) ...
    - 1/2*rhoLid.*uLid;
return
end



function [rho,us,vs] = CavityMacroBC(rho,us,vs, idx_xlo,idx_xhi,idx_ylo,idx_yhi,Ubc,Vbc)
us(idx_xlo) = 0; 
vs(idx_xlo) = 0;
us(idx_xhi) = 0;
vs(idx_xhi) = 0;
us(idx_ylo) = 0;
vs(idx_ylo) = 0;
us(idx_yhi) = Ubc(idx_yhi);
vs(idx_yhi) = 0;
return
end


function [ pEq ] = LBE_pEq(p,p0,u,v, qwgt,qex,qey,qc)
pEq = zeros(length(p),9);
vel2 = u.^2 + v.^2;
for i = 1:9
    eu = qc * (qex(i)*u + qey(i)*v);
    pEq(:,i) = qwgt(i)*(p + p0.*(3*eu + 9/2*(eu.^2) - 3/2*vel2));
    % pEq(:,i) = qwgt(i)*p.*(1 + 3*eu + 9/2*(eu.^2) - 3/2*vel2);
end
return
end

function [ pStream ] = LBE_pStream(p, nx,ny, qex,qey)
% pStream = zeros(size(p));
% for i = 1:9
    % pStream(:,i) = LBMStream2d(p(:,i), nx,ny, qex(i),qey(i));
% end
pStream = reshape(p,nx,ny,9);
i = 1; % fixed
%
i = 2; pStream(2:nx,:,i) = pStream(1:nx-1,:,i);
i = 3; pStream(:,2:ny,i) = pStream(:,1:ny-1,i);
i = 4; pStream(1:nx-1,:,i) = pStream(2:nx,:,i);
i = 5; pStream(:,1:ny-1,i) = pStream(:,2:ny,i);
%
i = 6; pStream(2:nx,2:ny,i) = pStream(1:nx-1,1:ny-1,i);
i = 7; pStream(1:nx-1,2:ny,i) = pStream(2:nx,1:ny-1,i);
i = 8; pStream(1:nx-1,1:ny-1,i) = pStream(2:nx,2:ny,i);
i = 9; pStream(2:nx,1:ny-1,i) = pStream(1:nx-1,2:ny,i);
pStream = reshape(pStream,nx*ny,9);
return
end

function [ pres,us,vs ] = LBEMacroVars(ps, pres0, qex,qey,qc)
pres = sum(ps,2);
us = qc * ps * qex' ./ pres0;
vs = qc * ps * qey' ./ pres0;
% us = qc * ps * qex' ./ pres;
% vs = qc * ps * qey' ./ pres;
return
end

function [ u,v ] = LBE_IBVelocity(u,v,ibu,ibv,ibfrac,ibidx)
if exist('ibidx','var')
    u(ibidx) = (1-ibfrac(ibidx)).*u(ibidx) + ibfrac(ibidx).*ibu(ibidx);
    v(ibidx) = (1-ibfrac(ibidx)).*v(ibidx) + ibfrac(ibidx).*ibv(ididx);
else
    u = (1-ibfrac).*u + ibfrac.*ibu;
    v = (1-ibfrac).*v + ibfrac.*ibv;
end
return
end

% if (0) % stiff
    % disp([int2str(step)]);
    % Kn = fcol;
    % Km = Kn;
    % iter_resid = zeros(1,9);
    % iter_conv = false;
    % for iter = 1:500
        % Kold = Km;
        
        % Km = 1/2 * (Kn + Kold);
        % relax = 1.0;
        % Km = relax*Km + (1-relax)*Kold;
        
        % % advect to f(m)
        % fsm = fs + omega*Km;
        % for i = 1:9
            % fm(:,i) = LBMStream2d(fsm(:,i), nx,ny, qex(i),qey(i));
        % end
        % fm = CavityMicroBC(fm,idx_xlo,idx_xhi,idx_ylo,idx_yhi,Ubc,Vbc,qc);
        % % macroscopic var.
        % [rhom,um,vm] = LBEMacroVars(fm, qex,qey,qc);
        % [rhom,um,vm] = CavityMacroBC(rhom,um,vm, idx_xlo,idx_xhi,idx_ylo,idx_yhi,Ubc,Vbc);
        
        % % new collision op.
        % feqm = LBMEqDistrib2D(rhom,um,vm,nx,ny, qwgt,qex,qey,qc,cs2);
        % Km = feqm - fm;
        
        % % check convergence
        % tol_rel = 1e-10;
        % for i = 1:9
            % iter_resid(i) = norm((Km(:,i)-Kold(:,i))./rhom, 2);
        % end
        % iter_conv = all(iter_resid < tol_rel);
        % if (iter_conv); break; end
    % end % end of relax iteration
    
    % if (iter_conv)
        % fstar = fs + omega*1/2*(Km+Km);
        % for i = 1:9
            % fs(:,i) = LBMStream2d(fstar(:,i), nx,ny, qex(i),qey(i));
        % end
    % else
        % disp(iter_resid);
        % error('Implicit Trapezoidal Iteration failed');
    % end
% else % trivial procedure
    % fstar = fs + omega*fcol;
    % % streaming
    % for i = 1:9
        % fs(:,i) = LBMStream2d(fstar(:,i), nx,ny, qex(i),qey(i));
    % end
    % % fstar = fs;
    % % for i = 1:9
        % % fs(:,i) = LBMStream2d(fstar(:,i), nx,ny, qex(i),qey(i));
        % % fs(:,i) = fs(:,i) + omega*fcol(:,i);
    % % end
% end


