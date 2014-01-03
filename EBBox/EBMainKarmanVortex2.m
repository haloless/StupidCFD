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

x_lo = 0;
x_hi = 4;
xlen = x_hi - x_lo;
y_lo = 0;
y_hi = 1;
ylen = y_hi - y_lo;
% refine = 10;
refine = 5;
nx = 32 * refine;
ny = 8 * refine;

rho = 1;
UIn = 1;
POut = 0;
L0 = ylen / 8;
Re = 100;
nu = UIn*L0 / Re;

object = 'cylinder';
switch object
    case {'square'}
        EBFracFunc = @(x,y) double(abs(x-1)<=(1/16) & abs(y-0.5)<=(1/16));
    case {'cylinder'}
        cylinder_D = L0;
        cylinder_R = cylinder_D/2;
        EBFracFunc = @(x,y) double((x-1).^2+(y-0.5).^2<=(cylinder_R)^2);
        EBDistFunc = @(x,y) sqrt((x-1).^2+(y-0.5).^2)-cylinder_R;
    otherwise
        error('Unknown object: %s', object);
end

dx = xlen / nx;
dy = ylen / ny;
dh = min([dx,dy]);

cellxs = linspace(x_lo-dx/2,x_hi+dx/2,nx+2);
cellys = linspace(y_lo-dy/2,y_hi+dy/2,ny+2);
edgexs = linspace(x_lo-dx,x_hi+dx,nx+3);
edgeys = linspace(y_lo-dy,y_hi+dy,ny+3);
[Xcell,Ycell] = ndgrid(cellxs,cellys);
%
Icell = 2:nx+1;
Jcell = 2:ny+1;
Iedge = 2:nx+2;
Jedge = 2:ny+2;

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


if exist('EBDistFunc','var')
    ain = 0.1*dh; aout = 0.8*dh;
    yc = ain / (ain+aout);
    
    ebls = EBDistFunc(Xcell,Ycell);
    ebvof = zeros(size(ebls));
    mask = ebls<=-ain; ebvof(mask) = 1;
    mask = ebls>-ain & ebls<=0; ebvof(mask) = yc + (yc-1)/ain*ebls(mask);
    mask = ebls>0 & ebls<=aout; ebvof(mask) = yc - yc/aout*ebls(mask);
    mask = ebls>aout; ebvof(mask) = 0;
    
    [Xedge,Yedge] = ndgrid(edgexs,cellys);
    ebls = EBDistFunc(Xedge,Yedge);
    ebvof_macx = zeros(size(ebls));
    mask = ebls<=-ain; ebvof_macx(mask) = 1;
    % mask = ebls>-ain & ebls<=0; ebvof_macx(mask) = yc + (yc-1)/ain*ebls(mask);
    mask = ebls>-ain & ebls<=0; ebvof_macx(mask) = 1;
    mask = ebls>0 & ebls<=aout; ebvof_macx(mask) = yc - yc/aout*ebls(mask);
    mask = ebls>aout; ebvof_macx(mask) = 0;
    
    [Xedge,Yedge] = ndgrid(cellxs,edgeys);
    ebls = EBDistFunc(Xedge,Yedge);
    ebvof_macy = zeros(size(ebls));
    mask = ebls<=-ain; ebvof_macy(mask) = 1;
    % mask = ebls>-ain & ebls<=0; ebvof_macy(mask) = yc + (yc-1)/ain*ebls(mask);
    mask = ebls>-ain & ebls<=0; ebvof_macy(mask) = 1;
    mask = ebls>0 & ebls<=aout; ebvof_macy(mask) = yc - yc/aout*ebls(mask);
    mask = ebls>aout; ebvof_macy(mask) = 0;
    
    clear ebls mask;
    clear Xedge Yedge;
else
    ebvof = EBFracFunc(Xcell,Ycell);
    %
    ebvof_macx = zeros(nx+3,ny+2);
    ebvof_macy = zeros(nx+2,ny+3);
    ebvof_macx(2:nx+2,2:ny+1) = 0.5*(ebvof(1:nx+1,2:ny+1) + ebvof(2:nx+2,2:ny+1));
    ebvof_macy(2:nx+1,2:ny+2) = 0.5*(ebvof(2:nx+1,1:ny+1) + ebvof(2:nx+1,2:ny+2));
end
% ebvof_macx(2:nx+2,2:ny+1) = 0.5*(ebvof(1:nx+1,2:ny+1) + ebvof(2:nx+2,2:ny+1));
% ebvof_macy(2:nx+1,2:ny+2) = 0.5*(ebvof(2:nx+1,1:ny+1) + ebvof(2:nx+1,2:ny+2));
%
eb_beta = 1 - ebvof;
eb_betax = 1 - ebvof_macx;
eb_betay = 1 - ebvof_macy;
% 
ebflag = (ebvof > 0);
ebflag_macx = (ebvof_macx > 0);
ebflag_macy = (ebvof_macy > 0);
%
ebflag_solid = (ebvof>0.999);

% initial condition
umac(Iedge,Jcell) = UIn;
% umac(Iedge,Jcell) = (1-ebvof_macx(Iedge,Jcell)) .* (UIn);

% time stepping
cfl = 0.25;
% dt_max = 1e-3;
dt_max = 5e-4;
dt = min([0.2*dh^2/nu, cfl*dh/UIn, dt_max]);

% build Poisson Op
disp('Building Poisson OP...');
tic;
% [LapOp rhs_corr] = PPELapOp(nx,ny,dx,dy);
[LapOp, rhs_corr] = EBPPELapOp(nx,ny,dx,dy, eb_betax(Iedge,Jcell),eb_betay(Icell,Jedge));
toc;
LapPerm = symamd(LapOp);
RLap = chol(LapOp(LapPerm,LapPerm)); RLapt = RLap';

% EBPPE iteration algo.
ebppe_A = zeros(nx,ny);
ebppe_Bx = eb_betax(Iedge,Jcell);
ebppe_By = eb_betay(Icell,Jedge);


% buffer for collecting fluid drag
% hydro_force = zeros(1,6);
hydro_force = [];


max_time = 50;
max_step = 800000;
% max_step = 10000;
time = 0;
step = 0;
while (time<max_time && step<max_step)
    time = time + dt;
    step = step + 1;
    
    [umac,vmac] = VelocityBC(umac,vmac,nx,ny);
    
    uold = umac;
    vold = vmac;
    pold = pres;
    
    % predictor
    [Cu,Cv,Du,Dv] = VelocityPredictor(umac,vmac,nx,ny,dx,dy,dt);
    Hu = Cu + Du;
    Hv = Cv + Dv;
    
    if (step==1)
        ustar = umac + dt*Hu;
        vstar = vmac + dt*Hv;
    else
        ustar = umac + dt/2*(3*Hu - Hu_old);
        vstar = vmac + dt/2*(3*Hv - Hv_old);
    end
    Hu_old = Hu;
    Hv_old = Hv;
    
    % EB direct forcing
    ustar(Iedge,Jcell) = (1-ebvof_macx(Iedge,Jcell)) .* ustar(Iedge,Jcell);
    vstar(Icell,Jedge) = (1-ebvof_macy(Icell,Jedge)) .* vstar(Icell,Jedge);
    
    % PPE
    % rhs = PPERhs(ustar,vstar,nx,ny,dx,dy,dt);
    rhs = EBPPERhs(ustar,vstar,nx,ny,dx,dy,dt, eb_betax,eb_betay);
    % rhs = rhs + rhs_corr;
    if (0)
        sol = rhs;
        sol(LapPerm) = RLap \ (RLapt \ rhs(LapPerm));
        pres(Icell,Jcell) = reshape(sol,nx,ny);
    else
        eps_rel = 1e-8;
        eps_abs = -1;
        max_sol_iter = 5000;
        sol = pres(Icell,Jcell);
        rhs = reshape(rhs,nx,ny);
        [ret,sol] = EBPPESolver_cg(nx,ny,dx,dy, ebppe_A,ebppe_Bx,ebppe_By, ...
        sol,rhs, eps_rel,eps_abs,max_sol_iter);
        if (ret~=0)
            error('EBPPE solver failure: errno=%d',ret);
        end
        pres(Icell,Jcell) = sol;
    end
    
    pres = PressureBC(pres,nx,ny);
    
    % corrector
    [ucorr,vcorr] = VelocityCorrector(pres,rho,nx,ny,dx,dy,dt);
    umac = ustar + dt*eb_betax.*ucorr;
    vmac = vstar + dt*eb_betay.*vcorr;
    [umac,vmac] = VelocityBC(umac,vmac,nx,ny);
    
    
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
            
            subnrow = 3;
            subncol = 2;
            
            subplot(subnrow,subncol,1);
            phi2 = pres;
            phi2(ebflag) = NaN;
            contourf(validxs,validys, phi2(2:nx+1,2:ny+1)',20);
            axis([x_lo x_hi y_lo y_hi]); 
            axis equal;
            title(['pressure ' prompt]);
            colorbar;
            
            subplot(subnrow,subncol,2);
            vel = sqrt(ucell.^2+vcell.^2);
            vel(ebflag) = NaN;
            contourf(validxs,validys, vel(2:nx+1,2:ny+1)', 20);
            axis([x_lo x_hi y_lo y_hi]); 
            axis equal;
            title(['velocity ']);
            colorbar;
            % hold on;
            % gap = 4;
            % quiver(cellxs(2:gap:nx+1),cellys(2:ny+1), ucell(2:gap:nx+1,2:ny+1)', vcell(2:gap:nx+1,2:ny+1)');
            % hold off;
            
            subplot(subnrow,subncol,3);
            psi = zeros(nx+2,ny+2);
            psi(2:nx+1,2:ny+1) = easy_streamfunc(validxs,validys,nx,ny,dx,dy, ...
                ucell(2:nx+1,2:ny+1), vcell(2:nx+1,2:ny+1));
            psi(ebflag) = NaN;
            contourf(validxs,validys, psi(2:nx+1,2:ny+1)', 20);
            axis([x_lo x_hi y_lo y_hi]);
            axis equal;
            title(['stream ']);
            colorbar;
            
            subplot(subnrow,subncol,4);
            vort = zeros(nx+2,ny+2);
            vort(2:nx+1,2:ny+1) = easy_vorticity(validxs,validys,nx,ny,dx,dy,ucell,vcell);
            % vort(2:nx+1,2:ny+1) = 1/(2*dx) * (vcell(3:nx+2,2:ny+1) - vcell(1:nx,2:ny+1)) ...
            % - 1/(2*dy) * (ucell(2:nx+1,3:ny+2) - ucell(2:nx+1,1:ny));
            vort(ebflag) = NaN;
            contourf(validxs,validys, vort(2:nx+1,2:ny+1)', 20);
            axis([x_lo x_hi y_lo y_hi]);
            axis equal;
            title(['vorticity ']);
            colorbar;
            
            subplot(subnrow,subncol,5);
            % drag and lift
            [diffu diffv] = VelocityDiffusion(uold,vold,nx,ny,dx,dy,dt);
            [gpx,gpy] = VelocityCorrector(pold,rho,nx,ny,dx,dy,dt);
            xdirf = rho * ebvof_macx .* (gpx + diffu);
            ydirf = rho * ebvof_macy .* (gpy + diffv);
            xdirf = 0.5 * (xdirf(1:nx+2,:)+xdirf(2:nx+3,:));
            ydirf = 0.5 * (ydirf(:,1:ny+2)+ydirf(:,2:ny+3));
            hydro_fx = dx*dy * sum(xdirf(ebflag));
            hydro_fy = dx*dy * sum(ydirf(ebflag));
            hydro_force(end+1,1:3) = [time, hydro_fx, hydro_fy];
            plot(hydro_force(:,1),hydro_force(:,2), hydro_force(:,1),hydro_force(:,3));
            legend('Drag (Fx)', 'Lift (Fy)');
            title(['hydrodynamics force ', num2str(hydro_fx)]);
            
            subplot(subnrow,subncol,6);
            plot(hydro_force(:,1),hydro_force(:,2)/(0.5*rho*UIn^2*cylinder_D), ...
                hydro_force(:,1),hydro_force(:,3)/(0.5*rho*UIn^2*cylinder_D));
            legend('C_drag (stress)', 'C_lift (stress)');
            title(['Cd=', num2str(hydro_fx/(0.5*rho*UIn^2*cylinder_D)), ...
            ';Cl=', num2str(hydro_fy/(0.5*rho*UIn^2*cylinder_D))]);
            
            % update head-up plots
            drawnow;
        end
    end
end






