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

% ## LBMMRT_CavityMain

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-16

clc; clear;

refine = 1;
Lx = 128 * refine;
Ly = 128 * refine;
nx = Lx + 2;
ny = Ly + 2;
hx = Lx / (nx-2);
hy = Ly / (ny-2);
hh = min([hx hy]);
qc = 1;
dt = hh / qc;
cs2 = 1/3 * qc^2;
cs = sqrt(cs2);

[X,Y] = ndgrid(linspace(-hx/2,Lx+hx/2,nx),linspace(-hy/2,Ly+hy/2,ny));

% determine parameters
Re = 400;
% U0 = 0.1;
omega = 1.6;

if exist('U0','var')
    nu = U0 * Lx / Re;
    tau = 1/cs2*nu/dt + 0.5;
    omega = 1/tau;
elseif exist('omega','var')
    tau = 1/omega;
    nu = (tau-0.5)*cs2*dt;
    U0 = Re*nu/Lx;
else
    error('Parameters must be given.');
end

% compute some indicators
Mach = U0 / cs;

% D2Q9
ex = [0, 1, 0, -1, 0, 1, -1, -1, 1];
ey = [0, 0, 1, 0, -1, 1, 1, -1, -1];
% MRT: s2,s3,s5,s7 to be chosen
% M: [1,rho; 2,e; 3,eps; 4,jx; 5,qx; 6,jy; 7,qy; 8,pxx; 9,pxy]
sRho = 1; sJx = 1; sJy = 1;
sE = 1.1; sEps = 1.1; sQx = 1.1; sQy = 1.1;
sPxx = 1/tau; sPxy = 1/tau;
S = [sRho, sE, sEps, sJx, sQx, sJy, sQy, sPxx, sPxy];
M = [1, 1, 1, 1, 1, 1, 1, 1, 1; ...
    -4,-1,-1,-1,-1, 2, 2, 2, 2; ...
    4,-2,-2,-2,-2, 1, 1, 1, 1; ...
    0, 1, 0,-1, 0, 1,-1,-1, 1; ...
    0,-2, 0, 2, 0, 1,-1,-1, 1; ...
    0, 0, 1, 0,-1, 1, 1,-1,-1; ...
    0, 0,-2, 0, 2, 1, 1,-1,-1; ...
    0, 1,-1, 1,-1, 0, 0, 0, 0; ...
    0, 0, 0, 0, 0, 1,-1, 1,-1];
InvMS = M \ diag(S);

% MRT coefficients, see [Lallemand & Luo, 2000]
alpha2 = -8; gamma2 = 18;
alpha3 = 4; gamma4 = -18;
c1 = -2;
gamma1 = 2/3; gamma3 = 2/3;
% mapping functions
meqRho = @(rho,jx,jy) rho;
meqE   = @(rho,jx,jy) 1/4*alpha2*rho + 1/6*gamma2*(jx.^2+jy.^2);
meqEps = @(rho,jx,jy) 1/4*alpha3*rho + 1/6*gamma4*(jx.^2+jy.^2);
meqJx  = @(rho,jx,jy) jx;
meqQx  = @(rho,jx,jy) 1/2*c1*jx;
meqJy  = @(rho,jx,jy) jy;
meqQy  = @(rho,jx,jy) 1/2*c1*jy;
meqPxx = @(rho,jx,jy) 3/2*gamma1 * (jx.^2-jy.^2);
meqPxy = @(rho,jx,jy) 3/2*gamma3 * (jx.*jy);

idx_xlo = sub2ind([nx,ny], 1*ones(1,ny), 1:ny);
idx_xhi = sub2ind([nx,ny], nx*ones(1,ny), 1:ny);
idx_ylo = sub2ind([nx,ny], 1:nx, 1*ones(1,nx));
idx_yhi = sub2ind([nx,ny], 1:nx, ny*ones(1,nx));

% storage
rho = ones(nx*ny,1);
u = zeros(nx*ny,1); u(idx_yhi) = U0;
v = zeros(nx*ny,1);
mEq = zeros(nx*ny,9);
f = zeros(nx*ny,9);
% initialize
mEq(:,1) = meqRho(rho,u,v);
mEq(:,2) = meqE(rho,u,v);
mEq(:,3) = meqEps(rho,u,v);
mEq(:,4) = meqJx(rho,u,v);
mEq(:,5) = meqQx(rho,u,v);
mEq(:,6) = meqJy(rho,u,v);
mEq(:,7) = meqQy(rho,u,v);
mEq(:,8) = meqPxx(rho,u,v);
mEq(:,9) = meqPxy(rho,u,v);
f = mEq / M';

max_step = 500000;
max_time = max_step * dt;
time = 0;
step = 0;
while (step<max_step && time<max_time)
    time = time + dt;
    step = step + 1;
    
    % collision
    jx = rho .* u;
    jy = rho .* v;
    mEq(:,1) = meqRho(rho,jx,jy);
    mEq(:,2) = meqE(rho,jx,jy);
    mEq(:,3) = meqEps(rho,jx,jy);
    mEq(:,4) = meqJx(rho,jx,jy);
    mEq(:,5) = meqQx(rho,jx,jy);
    mEq(:,6) = meqJy(rho,jx,jy);
    mEq(:,7) = meqQy(rho,jx,jy);
    mEq(:,8) = meqPxx(rho,jx,jy);
    mEq(:,9) = meqPxy(rho,jx,jy);
    
    m = f * M';
    fstar = reshape(f - (m-mEq)*InvMS', nx,ny,9);
    
    % streaming
    for i = 1:9
        f(:,i) = reshape(circshift(fstar(:,:,i), [ex(i),ey(i)]), nx*ny,1);
    end
    
    % microscopic BC
    % x-low wall, bounce back
    f(idx_xlo,2) = f(idx_xlo,4);
    f(idx_xlo,6) = f(idx_xlo,8);
    f(idx_xlo,9) = f(idx_xlo,7);
    % x-high wall, bounce back
    f(idx_xhi,4) = f(idx_xhi,2);
    f(idx_xhi,8) = f(idx_xhi,6);
    f(idx_xhi,7) = f(idx_xhi,9);
    % y-low wall, bounce back
    f(idx_ylo,3) = f(idx_ylo,5);
    f(idx_ylo,6) = f(idx_ylo,8);
    f(idx_ylo,7) = f(idx_ylo,9);
    % y-high wall, Zou-He
    rhoLid = f(idx_yhi,1) + f(idx_yhi,2) + f(idx_yhi,4) ...
        + 2 * (f(idx_yhi,3) + f(idx_yhi,6) + f(idx_yhi,7));
    % rhoLid = sum(f(idx_yhi,[1,2,4]),2) + 2*sum(f(idx_yhi,[3,6,7]),2);
    f(idx_yhi,5) = f(idx_yhi,3);
    f(idx_yhi,8) = f(idx_yhi,6) + 1/2*(f(idx_yhi,2)-f(idx_yhi,4)) - 1/2*rhoLid*U0;
    f(idx_yhi,9) = f(idx_yhi,7) + 1/2*(f(idx_yhi,4)-f(idx_yhi,2)) + 1/2*rhoLid*U0;
    
    
    % macroscopic
    uold = u;
    vold = v;
    
    rho = sum(f,2);
    u = f*ex' ./ rho;
    v = f*ey' ./ rho;
    
    % macroscopic BC
    u(idx_xlo) = 0; v(idx_xlo) = 0;
    u(idx_xhi) = 0; v(idx_xhi) = 0;
    u(idx_ylo) = 0; v(idx_ylo) = 0;
    u(idx_yhi) = U0; v(idx_yhi) = 0;
    
    % check steady
    tol_abs = 1e-7 * U0;
    resid = max([norm(u-uold,Inf), norm(v-vold,Inf)]);
    % resid = max([norm(u-uold,2), norm(v-vold,2);])
    is_conv = resid<=tol_abs;
    
    
    if (mod(step,100)==0 || is_conv)
        prompt = ['step=',int2str(step), ...
            ';Re=',num2str(Re), ';U0=',num2str(U0), ';Mach=',num2str(Mach), ...
            ';omega=',num2str(omega), ...
            ';|res/tol|=',num2str(resid/tol_abs)];
        disp(prompt);
        
        velx = reshape(u,nx,ny);
        vely = reshape(v,nx,ny);
        
        subplot(1,2,1);
        velmag = sqrt(velx.^2+vely.^2);
        contourf(X',Y', velmag' ./ U0, 32);
        colorbar;
        shading flat;
        axis equal; axis([0 Lx 0 Ly]);
        title([prompt]);
        
        subplot(1,2,2);
        psi = easy_streamfunc(X,Y,nx,ny,hx,hy,velx,vely);
        contourf(psi(2:nx-1,2:ny-1)',32);
        % vort = easy_vorticity(X,Y,nx-2,ny-2,hx,hy,u,v);
        % contourf(vort',32);
        colorbar;
        axis equal; axis([0 Lx 0 Ly]);
        
        drawnow;
    end
    
    if (is_conv); break; end
end % end of main loop

if (1 && Re==400)
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
    
    uprob = 1/U0 * 1/2*(velx(round(nx/2),:) + velx(round((nx+1)/2),:));
    vprob = 1/U0 * 1/2*(vely(:,round(ny/2)) + vely(:,round((ny+1)/2)));
    
    figure;
    plot(uprob,Y(1,:)./Ly,'-', GhiaU(:,1),GhiaU(:,2),'o'); 
    legend('LBM','Ghia'); title('U'); 
    axis equal; axis([-0.5 1 0 1]);
    
    figure;
    plot(X(:,1)./Lx,vprob,'-', GhiaV(:,1),GhiaV(:,2),'o'); 
    legend('LBM','Ghia'); title('V'); 
    axis equal; axis([0 1 -0.5 0.5])
end







