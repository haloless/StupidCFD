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
refine = 10;
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
%
pcorr = zeros(nx+2,ny+2);
% direct forcing
udirf = zeros(size(umac));
vdirf = zeros(size(vmac));


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
    mask = ebls>-ain & ebls<=0; ebvof_macx(mask) = yc + (yc-1)/ain*ebls(mask);
    mask = ebls>0 & ebls<=aout; ebvof_macx(mask) = yc - yc/aout*ebls(mask);
    mask = ebls>aout; ebvof_macx(mask) = 0;
    
    [Xedge,Yedge] = ndgrid(cellxs,edgeys);
    ebls = EBDistFunc(Xedge,Yedge);
    ebvof_macy = zeros(size(ebls));
    mask = ebls<=-ain; ebvof_macy(mask) = 1;
    mask = ebls>-ain & ebls<=0; ebvof_macy(mask) = yc + (yc-1)/ain*ebls(mask);
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
% 
ebflag = (ebvof > 0);
ebflag_macx = (ebvof_macx > 0);
ebflag_macy = (ebvof_macy > 0);
%
ebflag_solid = (ebvof>0.999);
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
umac(2:nx+2,2:ny+1) = UIn;
% time stepping
cfl = 0.25;
% dt_max = 1e-3;
dt_max = 5e-4;
dt = min([cfl*dh/UIn, dt_max]);

% build Poisson Op
disp('Building Poisson OP...');
tic;
[LapOp rhs_corr] = PPELapOp(nx,ny,dx,dy);
toc;
LapPerm = symamd(LapOp);
RLap = chol(LapOp(LapPerm,LapPerm)); RLapt = RLap';
% build Velocity OP
disp('Building Velocity OP...');
tic;
[LapU,LapV] = VelocityLapOp(nx,ny,dx,dy,dt);
toc;
LapUPerm = symamd(LapU);
RLapU = chol(LapU(LapUPerm,LapUPerm)); RLapUt = RLapU';
LapVPerm = symamd(LapV);
RLapV = chol(LapV(LapVPerm,LapVPerm)); RLapVt = RLapV';

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
    if (0)
        % impose EB
        umac(2:nx+2,2:ny+1) = (1-ebvof_macx(2:nx+2,2:ny+1)) .* umac(2:nx+2,2:ny+1);
        vmac(2:nx+1,2:ny+2) = (1-ebvof_macy(2:nx+1,2:ny+2)) .* vmac(2:nx+1,2:ny+2);
    end
    
    uold = umac;
    vold = vmac;
    pold = pres;
    
    % predictor
    [Hu,Hv,Du,Dv] = VelocityPredictor(umac,vmac,nx,ny,dx,dy,dt);
    [Gu,Gv] = VelocityCorrector(pres,rho,nx,ny,dx,dy,dt);
    
    if (1)
        % ustar(2:nx+2,2:ny+1) = umac(2:nx+2,2:ny+1) + dt*Hu(2:nx+2,2:ny+1) + dt*Du(2:nx+2,2:ny+1);
        % vstar(2:nx+1,2:ny+2) = vmac(2:nx+1,2:ny+2) + dt*Hv(2:nx+1,2:ny+2) + dt*Dv(2:nx+1,2:ny+2);
        ustar = umac + dt*(Hu+Du+Gu);
        vstar = vmac + dt*(Hv+Dv+Gv);
        [ustar,vstar] = VelocityBC(ustar,vstar,nx,ny);
    else
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
        
        ustar = ustar + dt*Gu;
        vstar = vstar + dt*Gv;
        
        [ustar,vstar] = VelocityBC(ustar,vstar,nx,ny);
        [urhs,vrhs] = VelocityLapRhs(ustar,vstar,ustar,vstar,nx,ny,dx,dy,dt);
        usol = urhs; usol(LapUPerm) = RLapU \ (RLapUt \ urhs(LapUPerm));
        vsol = vrhs; vsol(LapVPerm) = RLapV \ (RLapVt \ vrhs(LapVPerm));
        ustar(2:nx+2,2:ny+1) = reshape(usol,nx+1,ny);
        vstar(2:nx+1,2:ny+2) = reshape(vsol,nx,ny+1);
        [ustar,vstar] = VelocityBC(ustar,vstar,nx,ny);
    end
    
    if (1)
        % EB direct forcing
        udirf(2:nx+2,2:ny+1) = -ebvof_macx(2:nx+2,2:ny+1) .* ustar(2:nx+2,2:ny+1);
        vdirf(2:nx+1,2:ny+2) = -ebvof_macy(2:nx+1,2:ny+2) .* vstar(2:nx+1,2:ny+2);
        ustar(2:nx+2,2:ny+1) = ustar(2:nx+2,2:ny+1) + udirf(2:nx+2,2:ny+1);
        vstar(2:nx+1,2:ny+2) = vstar(2:nx+1,2:ny+2) + vdirf(2:nx+1,2:ny+2);
    end
    
    rhs = PPERhs(ustar,vstar,nx,ny,dx,dy,dt);
    rhs = rhs + rhs_corr;
    sol = rhs;
    sol(LapPerm) = RLap \ (RLapt \ rhs(LapPerm));
    
    pcorr(2:nx+1,2:ny+1) = reshape(sol,nx,ny);
    pcorr = PressureBC(pcorr,nx,ny);
    
    % corrector
    % [ucorr,vcorr] = VelocityCorrector(pres,rho,nx,ny,dx,dy,dt);
    [ucorr,vcorr] = VelocityCorrector(pcorr,rho,nx,ny,dx,dy,dt);
    umac(2:nx+2,2:ny+1) = ustar(2:nx+2,2:ny+1) + dt*ucorr(2:nx+2,2:ny+1);
    vmac(2:nx+1,2:ny+2) = vstar(2:nx+1,2:ny+2) + dt*vcorr(2:nx+1,2:ny+2);
    [umac,vmac] = VelocityBC(umac,vmac,nx,ny);
    %
    pres = pres + pcorr;
    pres = PressureBC(pres,nx,ny);

    
    % % EB direct forcing
    % udirf(2:nx+2,2:ny+1) = -ebvof_macx(2:nx+2,2:ny+1) .* umac(2:nx+2,2:ny+1);
    % vdirf(2:nx+1,2:ny+2) = -ebvof_macy(2:nx+1,2:ny+2) .* vmac(2:nx+1,2:ny+2);
    % % umac(2:nx+2,2:ny+1) = (1-ebvof_macx(2:nx+2,2:ny+1)) .* umac(2:nx+2,2:ny+1);
    % % vmac(2:nx+1,2:ny+2) = (1-ebvof_macy(2:nx+1,2:ny+2)) .* vmac(2:nx+1,2:ny+2);
    % umac(2:nx+2,2:ny+1) = umac(2:nx+2,2:ny+1) + udirf(2:nx+2,2:ny+1);
    % vmac(2:nx+1,2:ny+2) = vmac(2:nx+1,2:ny+2) + vdirf(2:nx+1,2:ny+2);
    
    if (mod(step,100)==0)
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
            % hydro_fx = -1/dt * rho*dx*dy * sum(udirf(ebflag_macx));
            % hydro_fy = -1/dt * rho*dx*dy * sum(vdirf(ebflag_macy));
            % xdirf = 0.5 * (udirf(1:nx+2,:) + udirf(2:nx+3,:));
            % ydirf = 0.5 * (vdirf(:,1:ny+2) + vdirf(:,2:ny+3));
            % hydro_fx = -1/dt * rho*dx*dy * sum(xdirf(ebflag));
            % hydro_fy = -1/dt * rho*dx*dy * sum(ydirf(ebflag));
            hydro_fxdf = -1/dt * rho*dx*dy * sum(udirf(ebflag_macx));
            hydro_fydf = -1/dt * rho*dx*dy * sum(vdirf(ebflag_macy));
            
            
            [diffu diffv] = VelocityDiffusion(uold,vold,nx,ny,dx,dy,dt);
            [gpx,gpy] = VelocityCorrector(pold,rho,nx,ny,dx,dy,dt);
            xdirf = rho * ebvof_macx .* (gpx + diffu);
            ydirf = rho * ebvof_macy .* (gpy + diffv);
            % hydro_fx1 = dx*dy * sum(xdirf(ebflag_macx));
            % hydro_fy1 = dx*dy * sum(ydirf(ebflag_macy));
            xdirf = 0.5 * (xdirf(1:nx+2,:)+xdirf(2:nx+3,:));
            ydirf = 0.5 * (ydirf(:,1:ny+2)+ydirf(:,2:ny+3));
            hydro_fx = dx*dy * sum(xdirf(ebflag));
            hydro_fy = dx*dy * sum(ydirf(ebflag));
            % hydro_fx2 = dx*dy * sum(xdirf(:));
            
            % hydro_force(end+1,1:3) = [time, hydro_fx, hydro_fy];
            % plot(hydro_force(:,1),hydro_force(:,2), hydro_force(:,1),hydro_force(:,3));
            % legend('Drag (Fx)', 'Lift (Fy)');
            % hydro_force(end+1,1:5) = [time,hydro_fx,hydro_fy,hydro_fx1,hydro_fx2];
            % plot(hydro_force(:,1),hydro_force(:,2), hydro_force(:,1),hydro_force(:,3), ...
                % hydro_force(:,1),hydro_force(:,4), hydro_force(:,1),hydro_force(:,5));
            % legend('Fx','Fy','Fx-MAC','Fx-total');
            % xlim([0 max_time]);
            hydro_force(end+1,1:5) = [time, hydro_fx, hydro_fy, hydro_fxdf,hydro_fydf];
            plot(hydro_force(:,1),hydro_force(:,2), hydro_force(:,1),hydro_force(:,3), ...
                hydro_force(:,1),hydro_force(:,4));
            legend('Drag (Fx)', 'Lift (Fy)', 'Drag (direct)');
            title(['hydrodynamics force ', num2str(hydro_fx)]);
            
            subplot(subnrow,subncol,6);
            plot(hydro_force(:,1),hydro_force(:,2)/(0.5*rho*UIn^2*cylinder_D), ...
                hydro_force(:,1),hydro_force(:,4)/(0.5*rho*UIn^2*cylinder_D), ...
                hydro_force(:,1),hydro_force(:,3)/(0.5*rho*UIn^2*cylinder_D));
            legend('Cd (stress)', 'Cd (direct)', 'Cl (stress)');
            title(['Cd=', num2str(hydro_fx/(0.5*rho*UIn^2*cylinder_D)), ...
            ';Cl=', num2str(hydro_fy/(0.5*rho*UIn^2*cylinder_D))]);
            
            % update head-up plots
            drawnow;
        end
    end
end






