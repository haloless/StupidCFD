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

function LBMSRT_CavityMain
clc;
clear all;

refine = 1;
Lx = 128 * refine;
Ly = 128 * refine;
nx = Lx + 1;
ny = Ly + 1;

hx = Lx / (nx-1);
hy = Ly / (ny-1);
hh = min([hx hy]);
qc = 1;
dt = hh / qc;
cs = 1/sqrt(3) * qc;
cs2 = 1/3 * qc^2;

% determine parameters
Re = 400;
% relaxation parameter, must < 2
omega = 1.7;
% ULid = 0.05;
% ULid = 0.1;

if (exist('omega','var'))
    tau = 1 / omega;
    nu = (tau-0.5) * cs2 * dt;
    ULid = Re*nu/Lx;
elseif (exist('ULid','var'))
    nu = ULid*Lx / Re;
    tau = 1/cs2*nu/dt + 0.5;
    omega = 1/tau;
else
    nu = 1e-3;
    ULid = Re * nu / Lx;
    tau = 1/cs2*nu/dt + 0.5;
    omega = 1/tau;
end

% compute some parameters
Re_h = ULid * hh / nu;
Cfl = dt*ULid / hh;
Dfl = dt*nu / hh^2;
Mach = ULid / cs;

% D2Q9
qwgt = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36];
qex  = [0, 1, 0, -1, 0, 1, -1, -1, 1];
qey  = [0, 0, 1, 0, -1, 1, 1, -1, -1];
qord = [1, 2, 3, 4, 5,  6, 7, 8,  9];
qopp = [1, 4, 5, 2, 3,  8, 9, 6,  7];


ILid = 2:(nx-1);
JLid = ny * ones(size(ILid));
IdxLid = sub2ind([nx,ny], ILid,JLid);

[X,Y] = ndgrid(linspace(0,Lx,nx),linspace(0,Ly,ny));

eb = ones(nx,ny);
eb(ILid,2:ny) = 0;
ebRegion = find(eb);

idx_xlo = sub2ind([nx,ny], 1*ones(1,ny), 1:ny);
idx_xhi = sub2ind([nx,ny], nx*ones(1,ny), 1:ny);
idx_ylo = sub2ind([nx,ny], 1:nx, 1*ones(1,nx));
idx_yhi = sub2ind([nx,ny], 1:nx, ny*ones(1,nx));

Ubc = zeros(nx*ny,1); Ubc(idx_yhi) = ULid;
Vbc = zeros(nx*ny,1); Vbc(idx_yhi) = 0;

% initial condition
fs = ones(nx*ny,1) * qwgt;
feq = zeros(nx*ny,9);
fstar = zeros(nx*ny,9);
% collision OP
fcol = zeros(nx*ny,9);

rho = ones(nx*ny,1);
us = zeros(nx*ny,1);
vs = zeros(nx*ny,1);

do_fracStep = 0;
if (do_fracStep)
    tauStar = 1;
    omegaStar = 1 / tauStar;
    nuStar = (tauStar-0.5) * cs2 * dt;
    
    FsOp = sparse(nx*ny,nx*ny);
    hx2 = hx^2; hy2 = hy^2;
    coefvec = dt*(nu-nuStar) * [-1/hy2, -1/hx2, 2/hx2+2/hy2, -1/hx2, -1/hy2] + [0 0 1 0 0];
    idx = 0;
    for j = 1:ny
        for i = 1:nx
            idx = idx + 1;
            
            if (j~=1 && j~=ny && i~=1 && i~=nx)
                FsOp(idx,[idx-nx,idx-1,idx,idx+1,idx+nx]) = coefvec;
            else
                FsOp(idx,idx) = 1;
            end
        end
    end
    
    eb = ones(nx,ny);
    eb(2:nx-1,2:ny-1) = 0;
    IdxWall = find(eb);
end
max_step = 1000000;
max_time = max_step*dt;
time = 0;
step = 0;
while (step<max_step && time<max_time)
    time = time + dt;
    step = step + 1;
    
    uold = us;
    vold = vs;
    
    % collision
    feq = LBMEqDistrib2D(rho,us,vs,nx,ny, qwgt,qex,qey, qc,cs2);
    fcol = feq - fs;
    
    if (do_fracStep)
        fstar = fs + omegaStar*fcol;
    else
        fstar = fs + omega*fcol;
    end
    % streaming
    for i = 1:9
        fs(:,i) = LBMStream2d(fstar(:,i), nx,ny, qex(i),qey(i));
    end
    
    % MICRO BC
    fs = CavityMicroBC(fs,idx_xlo,idx_xhi,idx_ylo,idx_yhi,Ubc,Vbc,qc);
    
    % MACRO
    [rho,us,vs] = LBEMacroVars(fs, qex,qey,qc);    
    % MACRO BC
    [rho,us,vs] = CavityMacroBC(rho,us,vs, idx_xlo,idx_xhi,idx_ylo,idx_yhi,Ubc,Vbc);
    
    if (do_fracStep)
        % rhs = us;
        % us = FsOp \ rhs;
        
        % rhs = vs;
        % vs = FsOp \ rhs;
        
        us = reshape(us,nx,ny);
        vs = reshape(vs,nx,ny);
        
        I = 2:nx-1;
        J = 2:ny-1;
        us(I,J) = us(I,J) + dt*(nu-nuStar) * ( ...
            1/hx^2*(us(I+1,J)-2*us(I,J)+us(I-1,J)) + 1/hy^2*(us(I,J+1)-2*us(I,J)+us(I,J-1)));
        vs(I,J) = vs(I,J) + dt*(nu-nuStar) * ( ...
            1/hx^2*(vs(I+1,J)-2*vs(I,J)+vs(I-1,J)) + 1/hy^2*(vs(I,J+1)-2*vs(I,J)+vs(I,J-1)));
        
        us = reshape(us,nx*ny,1);
        vs = reshape(vs,nx*ny,1);
        % [rho,us,vs] = CavityMacroBC(rho,us,vs, idx_xlo,idx_xhi,idx_ylo,idx_yhi,Ubc,Vbc);
    end
    
    % check convergence
    eps_abs = 1e-5 * ULid;
    res_abs = norm([us-uold; vs-vold],Inf);
    is_conv = res_abs<eps_abs;
    
    if (mod(step,100)==0 || is_conv)
        prompt = ['step=',int2str(step), ';time=',num2str(time), ...
        ';omega=',num2str(omega), ';tau=',num2str(tau), ';Mach=',num2str(Mach), ...
        ';Re=',num2str(Re), ';ULid=',num2str(ULid), ';nu=',num2str(nu), ...
        ';|res/tol|=',num2str(res_abs/eps_abs), ...
        ];
        disp(prompt);
        
        velx = reshape(us,nx,ny);
        vely = reshape(vs,nx,ny);
        subrow = 1;
        subcol = 2;
        
        % subplot(subrow,subcol,1);
        velmag = sqrt(velx.^2 + vely.^2);
        % velmag(ebRegion) = NaN;
        contourf(X',Y',velmag'./ULid, 32);
        colorbar; 
        shading flat;
        axis equal;
        title(prompt);
        
        % subplot(subrow,subcol,2);
        % psi = easy_streamfunc(X,Y,nx,ny,hx,hy,velx,vely);
        % psi(ebRegion) = NaN;
        % psi(Y>Ly*1/3) = NaN;
        % contourf(X',Y',psi',32);
        % colorbar;
        % axis equal;
        
        % [startx,starty] = meshgrid(round(nx/4),1:round(ny/32):ny);
        % % streamline(stream2(X',Y',velx',vely',startx,starty));
        % streamline(X',Y',velx',vely',startx,starty);
        % axis equal;
        % hold off;
        
        % stride = 4;
        % quiver(X(1:stride:nx,1:stride:ny)',Y(1:stride:nx,1:stride:ny)', ...
        % velx(1:stride:nx,1:stride:ny)',vely(1:stride:nx,1:stride:ny)', 2.0);
        % axis equal;
        % title([prompt]);
        
        % subplot(subrow,subcol,2);
        % pres = reshape(cs2*rho, nx,ny);
        % pres(ebRegion) = NaN;
        % contourf(X',Y',pres',16);
        % colorbar;
        % axis equal;
        
        drawnow;
        
        if (is_conv); break; end
    end
end % end of main loop

if (1)
    % velx = velx(validI,validJ);
    % vely = vely(validI,validJ);
    
    iprob = round(nx/2);
    jprob = round(ny/2);
    
    figure;
    plot(velx(iprob,:)./ULid,(1:ny)/ny);
    title('U'); axis equal; axis([-0.5 1 0 1])
    figure;
    plot((1:nx)/nx,vely(:,jprob)'./ULid);
    title('V'); axis equal; axis([0 1 -0.5 0.5]);
end
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

function [rho,us,vs] = LBEMacroVars(fs, qex,qey,qc)
rho = sum(fs,2);
us = qc * fs * qex' ./ rho;
vs = qc * fs * qey' ./ rho;
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

function [ K ] = CollisionOp(f, nx,ny, qwgt,qex,qey,qc,cs2, ...
idx_xlo,idx_xhi,idx_ylo,idx_yhi,Ubc,Vbc)

[rho,u,v] = LBEMacroVars(f,qex,qey,qc);
[rho,u,v] = CavityMacroBC(rho,u,v,idx_xlo,idx_xhi,idx_ylo,idx_yhi,Ubc,Vbc);

fEq = LBMEqDistrib2D(rho,u,v,nx,ny, qwgt,qex,qey,qc,cs2);
K = fEq - f;
return
end

function [ fn ] = MultiFracStep(f,N,omega,K, nx,ny, qwgt,qex,qey,qc,cs2, ...
idx_xlo,idx_xhi,idx_ylo,idx_yhi,Ubc,Vbc)
M = N * 2;
f2m = f;
for m = 1:N
    % 2*m+1
    for i = 1:9
        fcorr(:,i) = -LBMStream2d(f2m(:,i),nx,ny,-qex(i),-qey(i));
    end
    fcorr = fcorr + f2m + omega*K;
    %
    f2m1 = f2m + 1/M * fcorr;
    f2m1 = CavityMicroBC(f2m1,idx_xlo,idx_xhi,idx_ylo,idx_yhi,Ubc,Vbc,qc);
    K = CollisionOp(f2m1, nx,ny, qwgt,qex,qey,qc,cs2, ...
idx_xlo,idx_xhi,idx_ylo,idx_yhi,Ubc,Vbc);
    
    % 2*m+2
    for i = 1:9
        fcorr(:,i) = LBMStream2d(f2m1(:,i),nx,ny,qex(i),qey(i));
    end
    fcorr = fcorr - f2m1 + omega*K;
    %
    f2m2 = f2m1 + 1/M * fcorr;
    f2m2 = CavityMicroBC(f2m2,idx_xlo,idx_xhi,idx_ylo,idx_yhi,Ubc,Vbc,qc);
    K = CollisionOp(f2m2, nx,ny, qwgt,qex,qey,qc,cs2, ...
idx_xlo,idx_xhi,idx_ylo,idx_yhi,Ubc,Vbc);
    
    f2m = f2m2;
    % f2m = CavityMicroBC(f2m,idx_xlo,idx_xhi,idx_ylo,idx_yhi,Ubc,Vbc,qc);
end

fn = f2m;
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


