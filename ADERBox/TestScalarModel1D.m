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

% ## TestScalarModel1D

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-01-22

% clc
clear all

ader_order = 3;
% ader_order = 4;
% ader_order = 5;

ADERWENOGlobals1D;
ADERWENOInit1D(ader_order);

M = MDegree;
N = NPoint;
Ng = NGrow;

Nd = N * N;
Nvar = 1;

% 1D Gaussian Quadrature point
gpos = GausEta;
gwgt = GausWgt;

SpaceTimeDGGlobals1D;
SpaceTimeDGInit1D(Nvar,0);


% solve the model problem of LeVeque and Yee
% du/dt + a*du/dx = -nu*u*(u-1)*(u-1/2)
% we always have a = 1 and change the value of nu

xlo = 0.0;
xhi = 1.0;
xlen = xhi - xlo;

ncell = 100;
dx = xlen / ncell;

ng = Ng + 1;
nx = ng + ncell + ng;
lo = ng + 1;
hi = ng + ncell;
vrange = lo:hi;

cellxs = linspace(xlo-dx/2-dx*(ng-1), xhi+dx/2+dx*(ng-1), nx);
edgexs = linspace(xlo-dx*ng, xhi+dx*ng, nx+1);

us = zeros(nx,Nvar);
us(cellxs<=0.3) = 1;
us(cellxs>0.3) = 0;

% physical problem
% nu = 1.0;
% nu = 10.0;
% nu = 100.0;
nu = 1000.0;
alpha = 1.0;
% cfl = 0.75;
cfl = 0.5;
% NOTE too small DT could be diffusive for this problem
% dt = min([cfl*dx/alpha, 0.001]);
dt = cfl*dx/alpha;

% flux and source
f_func = @(u) (alpha.*u);
s_func = @(u) (-nu.*u.*(u-1).*(u-1/2));
fstar_func = @(u) (dt/dx*alpha.*u);
sstar_func = @(u) (-dt*nu.*u.*(u-1).*(u-1/2));
sprim_func = @(u) (-dt*nu*(3*u.^2 - 3*u + 1/2));

max_time = 0.3;
% max_step = 1;
max_step = 100000;
time = 0.0;
step = 0;

while (time<max_time && step<max_step)
    time = time + dt;
    step = step + 1;
    
    % apply BC
    for i = 1:ng
        us(lo-i) = us(lo);
        us(hi+i) = us(hi);
    end
    
    uold = us;
    fold = fstar_func(us);
    Sold = sstar_func(us);
    
    % cell WENO reconstruction
    wh0 = zeros(N*Nvar, nx);
    for i = lo-1:hi+1
        wh0(:,i) = ADERWENOReconstruct1D(us(i-M:i+M),dx);
    end
    
    if (0)
        plot(cellxs(vrange),wh0(1,vrange),'x', ...
        cellxs(vrange),wh0(2,vrange),'o', cellxs(vrange),wh0(3,vrange),'s');
    end
    
    qhat = zeros(Nd,nx);
    fhat = zeros(Nd,nx);
    shat = zeros(Nd,nx);
    % DG predictor
    dg_iter_max = 0;
    for i = lo-1:hi+1
        wh = wh0(:,i);
        
        w0 = zeros(Nd,1); % initial guess
        
        % first obtain initial guess based on cell average value
        % u^{n+1} = u^n - d(f*)/dxi + S*
        % drop flux term and implicitly treat source term we solve 
        % u^{n+1} - S^{*,n+1} - u^n = 0
        % alternatively, we use Crank-Nicolson for the source term
        u0 = uold(i);
        s0 = Sold(i);
        Fguess = @(u) (u - 0.5*sstar_func(u) - 0.5*s0 - u0);
        Jguess = @(u) (1 - 0.5*sprim_func(u));
        [ret,u1] = NewtonRaphsonSolve(Fguess,Jguess,u0);
        if (ret~=0)
            msg = ['Failed to obtain initial guess for cell i=',int2str(i)];
            error(msg);
        end
        % use for all components
        w0(:) = u1;
        % TODO a better way is to do this initial guess for each component
        
        % tol_rel = 1e-10;
        tol_rel = 0;
        tol_abs = 1e-14;
        max_iter = 200;
        
        F0w = F0Mat * wh;
        
        qh = w0;
        Jh = diag(sprim_func(qh));
        % MJh = MMat * Jh;
        fh = fstar_func(qh);
        sh = sstar_func(qh);
        res = K1Mat*qh + KxMat*fh - MMat*sh - F0w;
        resnorm = norm(res);
        tol = resnorm*tol_rel + tol_abs;
        
        niter = 0;
        while (resnorm>tol && niter<max_iter)
            niter = niter + 1;
            
            A = K1Mat - MMat*Jh;
            rhs = F0w - KxMat*fh + MMat*sh - MMat*Jh*qh;
            
            qn = A \ rhs;
            qh = qn;
            
            Jh = diag(sprim_func(qh));
            fh = fstar_func(qh);
            sh = sstar_func(qh);
            
            res = K1Mat*qh + KxMat*fh - MMat*sh - F0w;
            resnorm = norm(res);
        end
        dg_iter_max = max([dg_iter_max,niter]);
        if (0)
            disp(['ITER=',int2str(niter), ...
            ';TOL=',num2str(tol),';|RES|=',num2str(resnorm)]);
        end
        if (resnorm > tol)
            disp(['Failed to solve DG predictor for cell=',int2str(i)]);
            disp(['ITER=',int2str(niter), ...
            ';TOL=',num2str(tol),';|RES|=',num2str(resnorm)]);
            error('DG failure');
        end
        
        qhat(:,i) = qh;
        fhat(:,i) = f_func(qh);
        shat(:,i) = s_func(qh);
    end % end of DG predictor
    
    fbar = zeros(nx+1,1);
    for i = lo:hi+1
        il = i-1;
        ir = i;
        
        % Rusanov flux
        asig = abs(alpha);
        
        fpos = 0;
        fneg = 0;
        ql = qhat(:,il);
        fl = fhat(:,il);
        qr = qhat(:,ir);
        fr = fhat(:,ir);
        for p = 1:Nd
            [px,pt] = ind2sub([N,N],p);
            fpos = fpos + 0.5*(fl(p)+asig*ql(p)) * LagrPsi1(px) * gwgt(pt);
            fneg = fneg + 0.5*(fr(p)-asig*qr(p)) * LagrPsi0(px) * gwgt(pt);
        end
        
        fbar(i) = fpos + fneg;
    end
    
    sbar = zeros(nx,1);
    for i = lo:hi
        savg = 0;
        for p = 1:Nd
            [px,pt] = ind2sub([N,N],p);
            savg = savg + shat(p,i)*gwgt(px)*gwgt(pt);
        end
        sbar(i) = savg;
    end
    
    I = lo:hi;
    us(I) = uold(I) + dt/dx*(fbar(I)-fbar(I+1)) + dt*sbar(I);
    
    
    if (mod(step,5) == 0)
        prompt = ['ADER-WENO-DG', int2str(ader_order), ...
        ';step=',int2str(step),';time=',num2str(time), ...
        ';DGIterMax=',int2str(dg_iter_max)];
        disp(prompt);
        
        % plot(cellxs(I),us(I),'x-',cellxs(I),sbar(I),'o-');
        % legend('ubar','sbar');
        % title(prompt);
        
        subplot(1,2,1);
        plot(cellxs(I),us(I),'x-');
        legend('ubar');
        subplot(1,2,2);
        plot(cellxs(I),sbar(I),'x-');
        legend('sbar');
        
        drawnow;
    end
end






