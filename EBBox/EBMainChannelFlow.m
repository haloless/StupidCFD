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

% ## EBBox

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-24

clc;
clear all;

% global variables
EBGlobals;

rho = 1;
% nu = 1e-5;
nu = 1e-2;
% nu = 0.1;
UIn = 1;
POut = 0;


x_lo = 0;
x_hi = 6;
xlen = x_hi - x_lo;
y_lo = 0;
y_hi = 1;
ylen = y_hi - y_lo;

Re = UIn * ylen / nu;

refine = 4;
nx = 48 * refine;
ny = 8 * refine;

% square
% EBFracFunc = @(x,y) double(abs(x-1)<=(1/16) & abs(y-0.5)<=(1/16));

% cylinder
% EBFracFunc = @(x,y) double((x-1).^2+(y-0.5).^2<=(1/16)^2);
% EBDistFunc = @(x,y) sqrt((x-1).^2+(y-0.5).^2)-0.125;

EBFracFunc = @(x,y) zeros(size(x));

dx = xlen / nx;
dy = ylen / ny;
dh = min([dx,dy]);

cellxs = linspace(x_lo-dx/2,x_hi+dx/2,nx+2);
cellys = linspace(y_lo-dy/2,y_hi+dy/2,ny+2);
edgexs = linspace(x_lo-dx,x_hi+dx,nx+3);
edgeys = linspace(y_lo-dy,y_hi+dy,ny+3);
[Xcell,Ycell] = ndgrid(cellxs,cellys);

umac = zeros(nx+3,ny+2);
vmac = zeros(nx+2,ny+3);
pres = zeros(nx+2,ny+2);
% storage buffers
ustar = zeros(size(umac));
vstar = zeros(size(vmac));
Hu = zeros(size(umac));
Hv = zeros(size(vmac));
Hu_old = Hu;
Hv_old = Hv;
Du = zeros(size(umac));
Dv = zeros(size(vmac));


ebvof = EBFracFunc(Xcell,Ycell);
% ebls = EBDistFunc(Xcell,Ycell);
ebflag = (ebvof > 0);
ebvof_macx = zeros(nx+3,ny+2);
ebvof_macy = zeros(nx+2,ny+3);
ebvof_macx(2:nx+2,2:ny+1) = 0.5*(ebvof(1:nx+1,2:ny+1) + ebvof(2:nx+2,2:ny+1));
ebvof_macy(2:nx+1,2:ny+2) = 0.5*(ebvof(2:nx+1,1:ny+1) + ebvof(2:nx+1,2:ny+2));
if (0)
    figure;
    subplot(2,1,1);
    contourf(Xcell',Ycell',ebvof');
    title('EB fraction');
    subplot(2,1,2);
    contourf(Xcell',Ycell',ebls');
    title('EB levelset');
    return
end

% initial condition
% umac(:,:) = UIn;
umac(2:nx+2,2:ny+1) = UIn;

cfl = 0.2;
dt_max = 5e-3;
dt = min([cfl*dh/UIn, dt_max]);

% build Poisson OP
disp('Building Poisson OP...');
tic
[LapOp,rhs_corr] = PPELapOp(nx,ny,dx,dy);
toc
LapPerm = symamd(LapOp);
RLap = chol(LapOp(LapPerm,LapPerm));
RLapt = RLap';
% build Velocity OP
disp('Building Velocity OP...');
tic;
[LapU,LapV] = VelocityLapOp(nx,ny,dx,dy,dt);
toc;
LapUPerm = symamd(LapU);
RLapU = chol(LapU(LapUPerm,LapUPerm)); RLapUt = RLapU';
LapVPerm = symamd(LapV);
RLapV = chol(LapV(LapVPerm,LapVPerm)); RLapVt = RLapV';

%
ppe_A = zeros(nx,ny);
ppe_Bx = ones(nx+1,ny);
ppe_By = ones(nx,ny+1);

max_time = 8;
max_step = 200000;
% max_step = 1;
time = 0;
step = 0;
while (time<max_time && step<max_step)
    time = time + dt;
    step = step + 1;
    
    [umac,vmac] = VelocityBC(umac,vmac,nx,ny);
    % impose EB
    umac(2:nx+2,2:ny+1) = (1-ebvof_macx(2:nx+2,2:ny+1)) .* umac(2:nx+2,2:ny+1);
    vmac(2:nx+1,2:ny+2) = (1-ebvof_macy(2:nx+1,2:ny+2)) .* vmac(2:nx+1,2:ny+2);
    
    uold = umac;
    vold = vmac;
    
    % predictor
    [Hu,Hv,Du,Dv] = VelocityPredictor(umac,vmac,nx,ny,dx,dy,dt);
    
    % ustar(2:nx+2,2:ny+1) = umac(2:nx+2,2:ny+1) + dt*Hu(2:nx+2,2:ny+1);
    % vstar(2:nx+1,2:ny+2) = vmac(2:nx+1,2:ny+2) + dt*Hv(2:nx+1,2:ny+2);
    
    if (step == 1)
        ustar = umac + dt*Hu + dt/2*Du;
        vstar = vmac + dt*Hv + dt/2*Dv;
    else
        ustar = umac + dt/2*(3*Hu-Hu_old) + dt/2*Du;
        vstar = vmac + dt/2*(3*Hv-Hv_old) + dt/2*Dv;
    end
    %
    Hu_old = Hu;
    Hv_old = Hv;
    
    %
    [ustar,vstar] = VelocityBC(ustar,vstar,nx,ny);
    % urhs = reshape(ustar(2:nx+2,2:ny+1),(nx+1)*ny,1) + corrURhs;
    % vrhs = reshape(vstar(2:nx+1,2:ny+2),nx*(ny+1),1) + corrVRhs;
    % usol = LapU \ urhs;
    % vsol = LapV \ vrhs;
    [urhs,vrhs] = VelocityLapRhs(ustar,vstar,ustar,vstar,nx,ny,dx,dy,dt);
    usol = urhs; usol(LapUPerm) = RLapU \ (RLapUt \ urhs(LapUPerm));
    vsol = vrhs; vsol(LapVPerm) = RLapV \ (RLapVt \ vrhs(LapVPerm));
    ustar(2:nx+2,2:ny+1) = reshape(usol,nx+1,ny);
    vstar(2:nx+1,2:ny+2) = reshape(vsol,nx,ny+1);
    [ustar,vstar] = VelocityBC(ustar,vstar,nx,ny);
    
    
    rhs = PPERhs(ustar,vstar,nx,ny,dx,dy,dt);
    % rhs = rhs + rhs_corr;
    if (0)
        sol = rhs;
        sol(LapPerm) = RLap \ (RLapt \ rhs(LapPerm));
        pres(2:nx+1,2:ny+1) = reshape(sol,nx,ny);
    else
        sol = zeros(nx,ny);
        rhs = reshape(rhs,nx,ny);
        if (0)
            eps_rel = 1e-9;
            eps_abs = -1;
            max_sol_iter = 8000;
            [ret,sol] = EBPPESolver_cg(nx,ny,dx,dy, ppe_A,ppe_Bx,ppe_By, ...
                sol,rhs, eps_rel,eps_abs,max_sol_iter);
        else
            % params.precond = 'None';
            params.precond = 'jacobi';
            params.tol = 1e-8;
            %
            % method = 'gmres';
            % params.restart = 100;
            % params.maxit = 500;
            method = 'pcg';
            params.maxit = 8000;
            % method = 'bicgstab';
            % params.maxit = 8000;
            [ret,sol] = EBPPESolver(nx,ny,dx,dy, ppe_A,ppe_Bx,ppe_By, ...
                sol,rhs, method,params);
        end
        if (ret~=0)
            error('EBPPE solver failure: errno=%d',ret);
        end
        
        pres(2:nx+1,2:ny+1) = sol;
    end
    
    pres = PressureBC(pres,nx,ny);
    
    % corrector
    [ucorr,vcorr] = VelocityCorrector(pres,rho,nx,ny,dx,dy,dt);
    umac(2:nx+2,2:ny+1) = ustar(2:nx+2,2:ny+1) + dt*ucorr(2:nx+2,2:ny+1);
    vmac(2:nx+1,2:ny+2) = vstar(2:nx+1,2:ny+2) + dt*vcorr(2:nx+1,2:ny+2);
    
    [umac,vmac] = VelocityBC(umac,vmac,nx,ny);
    
    % EB direct forcing
    umac(2:nx+2,2:ny+1) = (1-ebvof_macx(2:nx+2,2:ny+1)) .* umac(2:nx+2,2:ny+1);
    vmac(2:nx+1,2:ny+2) = (1-ebvof_macy(2:nx+1,2:ny+2)) .* vmac(2:nx+1,2:ny+2);
    
    if (mod(step,10)==0)
        ucell = 0.5 * (umac(1:nx+2,:) + umac(2:nx+3,:));
        vcell = 0.5 * (vmac(:,1:ny+2) + vmac(:,2:ny+3));
        
        dumac = umac - uold;
        dvmac = vmac - vold;
        derr = max([norm(dumac(:),Inf), norm(dvmac(:),Inf)]);
        
        prompt = ['step=',int2str(step), '; time=',num2str(time), '; err=',num2str(derr)];
        disp(prompt);
        
        if (1)
            validxs = cellxs(2:nx+1);
            validys = cellys(2:ny+1);
            
            subplot(2,2,1);
            phi2 = pres;
            phi2(ebflag) = NaN;
            contourf(validxs,validys, phi2(2:nx+1,2:ny+1)',20);
            axis([x_lo x_hi y_lo y_hi]); 
            % axis equal;
            title(['pressure ' prompt]);
            colorbar;
            
            subplot(2,2,2);
            vel = sqrt(ucell.^2+vcell.^2);
            vel(ebflag) = NaN;
            contourf(validxs,validys, vel(2:nx+1,2:ny+1)', 20);
            axis([x_lo x_hi y_lo y_hi]); 
            % axis equal;
            title(['velocity ']);
            colorbar;
            % hold on;
            % gap = 4;
            % quiver(cellxs(2:gap:nx+1),cellys(2:ny+1), ucell(2:gap:nx+1,2:ny+1)', vcell(2:gap:nx+1,2:ny+1)');
            % hold off;
            
            subplot(2,2,3);
            psi = zeros(nx+2,ny+2);
            psi(2:nx+1,2:ny+1) = easy_streamfunc(validxs,validys,nx,ny,dx,dy, ...
                ucell(2:nx+1,2:ny+1), vcell(2:nx+1,2:ny+1));
            psi(ebflag) = NaN;
            contourf(validxs,validys, psi(2:nx+1,2:ny+1)', 20);
            axis([x_lo x_hi y_lo y_hi]);
            % axis equal;
            title(['stream ']);
            colorbar;
            
            subplot(2,2,4);
            vort = zeros(nx+2,ny+2);
            vort(2:nx+1,2:ny+1) = easy_vorticity(validxs,validys,nx,ny,dx,dy,ucell,vcell);
            % vort(2:nx+1,2:ny+1) = 1/(2*dx) * (vcell(3:nx+2,2:ny+1) - vcell(1:nx,2:ny+1)) ...
            % - 1/(2*dy) * (ucell(2:nx+1,3:ny+2) - ucell(2:nx+1,1:ny));
            vort(ebflag) = NaN;
            contourf(validxs,validys, vort(2:nx+1,2:ny+1)', 20);
            axis([x_lo x_hi y_lo y_hi]);
            % axis equal;
            title(['vorticity ']);
            colorbar;
            
            drawnow;
        end
    end
end

if (1)
    ys = cellys(2:ny+1);
    % exact value
    uexact = 6*UIn/ylen^2 * ys .* (ylen-ys);
    pdrop = 12*UIn/ylen^2 * (1/Re);
    
    % sample positions
    ia = 2;
    ua = umac(ia,2:ny+1);
    ib = nx + 2;
    ub = umac(ib,2:ny+1);
    ic = round((x_lo+xlen*0.5)/dx);
    uc = umac(ic,2:ny+1);
    id = round((x_lo+xlen*0.25)/dx);
    ud = umac(id,2:ny+1);
    
    % figure;
    % plot(ys,uexact,'-', ys,ua,'x', ys,ub,'o', ys,uc,'+');
    % legend('exact', '0', '1', '1/2');
    
    figure;
    plot(uexact,ys,'-', ua,ys,'x', ub,ys,'o', uc,ys,'+', ud,ys,'.');
    legend('exact', '0', '1', '1/2', '1/4');
end




