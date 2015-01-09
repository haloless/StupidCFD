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

x_lo = -2.0;
x_hi = 2.0;
xlen = x_hi - x_lo;
y_lo = -2.0;
y_hi = 2.0;
ylen = y_hi - y_lo;
% nx = 16;
% ny = 16;
% nx = 32;
% ny = 32;
% nx = 48;
% ny = 48;
nx = 64;
ny = 64;
% nx = 80;
% ny = 80;
% nx = 25;
% ny = 25;
% nx = 50;
% ny = 50;
% nx = 100;
% ny = 100;
% nx = 128;
% ny = 128;

rho = 1;
UIn = 0;
POut = 0;

% Couette flow
R0 = 0.5;
Omega0 = 0.1;
R1 = 1.5;
Omega1 = 0.0;

% determine viscosity
% Re = 40;
% Re = 10;
Re = 1;
% Re = 200;
nu = max(R0*abs(Omega0),R1*abs(Omega1)) * R0 / Re

nu_s = nu * 1e2
nu_s = nu


% setup grid
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

ebls = ProbDistFunc(Xcell,Ycell);
%
ebls_umac = zeros(nx+3,ny+2);
for i = 1:nx+3
for j = 1:ny+2
    ebls_umac(i,j) = ProbDistFunc(edgexs(i),cellys(j));
end
end
%
ebls_vmac = zeros(nx+2,ny+3);
for i = 1:nx+2
for j = 1:ny+3
    ebls_vmac(i,j) = ProbDistFunc(cellxs(i),edgeys(j));
end
end
%
ebvof = zeros(nx+2,ny+2);
for i = 1:nx+2
for j = 1:ny+2
    ebvof(i,j) = EBSampleCellFrac(Xcell(i,j),Ycell(i,j),dx,dy);
end
end

thickness = 1.0 * dh;
% thickness = 1.5 * dh;
%
ebvof_macx = zeros(nx+3,ny+2);
eb_umac = zeros(nx+3,ny+2);
for i = 1:nx+3
for j = 1:ny+2
    r = sqrt(edgexs(i)^2 + cellys(j)^2);
    if (R0<r && r<R1)
        ebvof_macx(i,j) = 0.0;
    else
        ebvof_macx(i,j) = 1.0;
    end
    if r < 0.5*(R0+R1)
        eb_umac(i,j) = -Omega0 * cellys(j);
    end
    
    % ebvof_macx(i,j) = EBSampleCellFrac(edgexs(i),cellys(j),dx,dy);
    
    d = ebls_umac(i,j);
    if d < -thickness
        ebvof_macx(i,j) = 1.0;
    elseif d > thickness
        ebvof_macx(i,j) = 0.0;
    else
        ebvof_macx(i,j) = 0.5 * (1.0 - d/thickness - 1.0/pi*sin(pi*d/thickness));
    end
end
end
%
ebvof_macy = zeros(nx+2,ny+3);
eb_vmac = zeros(nx+2,ny+3);
for i = 1:nx+2
for j = 1:ny+3
    r = sqrt(cellxs(i)^2 + edgeys(j)^2);
    if (R0<r && r<R1)
        ebvof_macy(i,j) = 0.0;
    else
        ebvof_macy(i,j) = 1.0;
    end
    if r < 0.5*(R0+R1)
        eb_vmac(i,j) = Omega0 * cellxs(i);
    end
    
    % ebvof_macy(i,j) = EBSampleCellFrac(cellxs(i),edgeys(j),dx,dy);
    
    d = ebls_vmac(i,j);
    if d < -thickness
        ebvof_macy(i,j) = 1.0;
    elseif d > thickness
        ebvof_macy(i,j) = 0.0;
    else
        ebvof_macy(i,j) = 0.5 * (1.0 - d/thickness - 1.0/pi*sin(pi*d/thickness));
    end
end
end

% 
ebflag = (ebvof > 0);
ebflag_macx = (ebvof_macx > 0);
ebflag_macy = (ebvof_macy > 0);
%
ebflag_solid = (ebvof>0.999);

if (0)
    figure;
    subplot(2,2,1);
    contourf(Xcell',Ycell',ebvof');
    title('EB fraction');
    subplot(2,2,2);
    contourf(Xcell',Ycell',ebls');
    title('EB levelset');
    subplot(2,2,3);
    % eb_ucell = eb_umac .* ebvof_macx;
    % eb_ucell = eb_umac;
    % eb_ucell = 0.5*(eb_ucell(1:nx+2,:)+eb_ucell(2:nx+3,:));
    % eb_vcell = eb_vmac .* ebvof_macy;
    % eb_vcell = eb_vmac;
    % eb_vcell = 0.5*(eb_vcell(:,1:ny+2)+eb_vcell(:,2:ny+3));
    % quiver(Xcell',Ycell', eb_ucell', eb_vcell');
    
    eb_unode = 0.5 * (eb_umac(2:nx+2,1:ny+1) + eb_umac(2:nx+2,2:ny+2));
    eb_vnode = 0.5 * (eb_vmac(1:nx+1,2:ny+2) + eb_vmac(2:nx+2,2:ny+2));
    quiver(edgexs(2:nx+2),edgeys(2:ny+2), eb_unode',eb_vnode');
    
    title('EB velocity');
    return
end

% initial condition
if (1)
    aa = (Omega1*R1^2 - Omega0*R0^2) / (R1^2 - R0^2);
    bb = (Omega0-Omega1) * R1^2 * R0^2 / (R1^2 - R0^2);
    %
    uana = zeros(nx+3,ny+2);
    for i = 1:nx+3
    for j = 1:ny+2
        r = sqrt(edgexs(i)^2 + cellys(j)^2);
        if (r <= R0)
            uana(i,j) = -Omega0 * cellys(j);
        elseif R0<r && r<R1
            uana(i,j) = -(aa*r + bb/r) * (cellys(j)/r);
        end
    end
    end
    %
    vana = zeros(nx+2,ny+2);
    for i = 1:nx+2
    for j = 1:ny+3
        r = sqrt(cellxs(i)^2 + edgeys(j)^2);
        if (r <= R0)
            vana(i,j) = Omega0 * cellxs(i);
        elseif R0<r && r<R1
            vana(i,j) = (aa*r + bb/r) * (cellxs(i)/r);
        end
    end
    end
    %
    pp = zeros(nx+2,ny+2);
    kk = Omega0 * R0^2 / (R1^2-R0^2);
    r = R0;
    pp0 = kk^2 * (0.5*r^2 - 0.5*R1^4/r^2 - R1^2*log(r^2));
    r = R1;
    pp1 = kk^2 * (0.5*r^2 - 0.5*R1^4/r^2 - R1^2*log(r^2));
    for i = 1:nx+2
    for j = 1:ny+2
        r = sqrt(cellxs(i)^2 + cellys(j)^2);
        if r <= R0
            pp(i,j) = pp0 - pp1;
        elseif R0<r && r<R1
            pp(i,j) = kk^2 * (0.5*r^2 - 0.5*R1^4/r^2 - R1^2*log(r^2)) - pp1;
        else
            pp(i,j) = 0;
        end
    end
    end
end
umac = (1.0-ebvof_macx).*umac + ebvof_macx.*eb_umac;
vmac = (1.0-ebvof_macy).*vmac + ebvof_macy.*eb_vmac;
% time stepping
cfl = 0.25;
% dt_max = 1e-3;
% dt_max = 5e-3;
dt_max = 10e-3;
% dt_max = 2.5e-4;
% dt_max = 1.0e-4;
% dt = min([cfl*dh/UIn, dt_max]);
dt = dt_max;

% build Poisson Op
disp('Building Poisson OP...');
tic;
[LapOp rhs_corr] = PPELapOp(nx,ny,dx,dy);
% set reference
LapOp(1,:) = 0;
LapOp(:,1) = 0;
LapOp(1,1) = 1;
toc;
LapPerm = symamd(LapOp);
RLap = chol(LapOp(LapPerm,LapPerm)); RLapt = RLap';

% build Velocity OP
disp('Building Velocity OP...');
tic;
LapVel = PMVelocityLapOp(ebvof, nx,ny,dx,dy,dt);
toc;
LapVelPerm = symamd(LapVel);
RLapVel = chol(LapVel(LapVelPerm,LapVelPerm));
RLapVelt = RLapVel';
%
udof = (1:(nx-1)*ny);
vdof = (1:nx*(ny-1)) + (nx-1)*ny;



max_time = 5000;
% max_step = 800000;
% max_step = 10;
max_step = 10000;
% max_step = 1;
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
    [Hu,Hv,Du,Dv] = VelocityConvection(umac,vmac,nx,ny,dx,dy,dt);
    % [Gu,Gv] = VelocityCorrector(pres,rho,nx,ny,dx,dy,dt);
    
    
    
    % ustar = umac + dt*(Hu+Du);
    % vstar = vmac + dt*(Hv+Dv);
    
    ustar = umac + dt*(Hu);
    vstar = vmac + dt*(Hv);
    
    
    urhs = reshape(ustar(3:nx+1,2:ny+1), (nx-1)*ny,1);
    vrhs = reshape(vstar(2:nx+1,3:ny+1), nx*(ny-1),1);
    velrhs = [urhs; vrhs];
    %
    velsol = zeros((nx-1)*ny+nx*(ny-1),1);
    velsol(LapVelPerm) = RLapVel \ (RLapVelt \ velrhs(LapVelPerm));
    %
    ustar(3:nx+1,2:ny+1) = reshape(velsol(udof), nx-1,ny);
    vstar(2:nx+1,3:ny+1) = reshape(velsol(vdof), nx,ny-1);
    
    [ustar,vstar] = VelocityBC(ustar,vstar,nx,ny);
    
    
    rhs = PPERhs(ustar,vstar,nx,ny,dx,dy,dt);
    % rhs0 = rhs;
    
    rhs(1) = 0;
    sol = rhs;
    sol(LapPerm) = RLap \ (RLapt \ rhs(LapPerm));
    
    pcorr(2:nx+1,2:ny+1) = reshape(sol,nx,ny);
    % pcorr(2:nx+1,2:ny+1) = reshape(sol,nx,ny) + dt*nucell(2:nx+1,2:ny+1).*reshape(rhs0,nx,ny);
    pcorr = PressureBC(pcorr,nx,ny);
    
    % corrector
    [ucorr,vcorr] = VelocityCorrector(pcorr,rho,nx,ny,dx,dy,dt);
    % umac(2:nx+2,2:ny+1) = ustar(2:nx+2,2:ny+1) + dt*LapUmac.*ucorr(2:nx+2,2:ny+1);
    % vmac(2:nx+1,2:ny+2) = vstar(2:nx+1,2:ny+2) + dt*LapVmac.*vcorr(2:nx+1,2:ny+2);
    umac(2:nx+2,2:ny+1) = ustar(2:nx+2,2:ny+1) + dt*ucorr(2:nx+2,2:ny+1);
    vmac(2:nx+1,2:ny+2) = vstar(2:nx+1,2:ny+2) + dt*vcorr(2:nx+1,2:ny+2);
    [umac,vmac] = VelocityBC(umac,vmac,nx,ny);
    %
    % pres = pres + pcorr;
    pres = pcorr;
    pres = PressureBC(pres,nx,ny);
    
    umac = (1.0-ebvof_macx).*umac + ebvof_macx.*eb_umac;
    vmac = (1.0-ebvof_macy).*vmac + ebvof_macy.*eb_vmac;
    
    if (mod(step,50)==0)
        unode = 0.5 * (umac(2:nx+2,1:ny+1) + umac(2:nx+2,2:ny+2));
        vnode = 0.5 * (vmac(1:nx+1,2:ny+2) + vmac(2:nx+2,2:ny+2));
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
            
            subnrow = 1;
            subncol = 2;
            
            subplot(subnrow,subncol,1);
            phi2 = pres;
            % phi2(ebflag) = NaN;
            % phi2(ebls<=0) = NaN;
            % contourf(validxs,validys, phi2(2:nx+1,2:ny+1)');
            imagesc(validxs,validys, phi2(2:nx+1,2:ny+1)');
            set(gca,'YDir','normal');
            axis([x_lo x_hi y_lo y_hi]); 
            axis equal;
            title(['pressure ' prompt]);
            colorbar;
            hold on;
            contour(validxs,validys, ebls(2:nx+1,2:ny+1)', [0,0],'LineColor','b','LineWidth',2);
            hold off;
            
            subplot(subnrow,subncol,2);
            vel = sqrt(ucell.^2+vcell.^2);
            % vel(ebflag) = NaN;
            contourf(validxs,validys, vel(2:nx+1,2:ny+1)');
            % imagesc(validxs,validys, vel(2:nx+1,2:ny+1)');
            axis([x_lo x_hi y_lo y_hi]); 
            axis equal;
            title(['velocity ']);
            colorbar;
            hold on;
            contour(validxs,validys, ebls(2:nx+1,2:ny+1)', ...
            [0.0,0.0],'LineColor','w','LineWidth',2);
            hold off;
            hold on;
            gap = 1;
            % quiver(cellxs(2:gap:nx+1),cellys(2:ny+1), ...
            % ucell(2:gap:nx+1,2:ny+1)', vcell(2:gap:nx+1,2:ny+1)');
            quiver(edgexs(2:nx+2),edgeys(2:ny+2), unode', vnode');
            hold off;
            
            % subplot(subnrow,subncol,3);
            % psi = zeros(nx+2,ny+2);
            % psi(2:nx+1,2:ny+1) = easy_streamfunc(validxs,validys,nx,ny,dx,dy, ...
                % ucell(2:nx+1,2:ny+1), vcell(2:nx+1,2:ny+1));
            % psi(ebflag) = NaN;
            % contourf(validxs,validys, psi(2:nx+1,2:ny+1)', 20);
            % axis([x_lo x_hi y_lo y_hi]);
            % axis equal;
            % title(['stream ']);
            % colorbar;
            
            % subplot(subnrow,subncol,4);
            % vort = zeros(nx+2,ny+2);
            % vort(2:nx+1,2:ny+1) = easy_vorticity(validxs,validys,nx,ny,dx,dy,ucell,vcell);
            % % vort(2:nx+1,2:ny+1) = 1/(2*dx) * (vcell(3:nx+2,2:ny+1) - vcell(1:nx,2:ny+1)) ...
            % % - 1/(2*dy) * (ucell(2:nx+1,3:ny+2) - ucell(2:nx+1,1:ny));
            % vort(ebflag) = NaN;
            % contourf(validxs,validys, vort(2:nx+1,2:ny+1)', 20);
            % axis([x_lo x_hi y_lo y_hi]);
            % axis equal;
            % title(['vorticity ']);
            % colorbar;
            
            % subplot(subnrow,subncol,5);
            % % hydro_fx = -1/dt * rho*dx*dy * sum(udirf(ebflag_macx));
            % % hydro_fy = -1/dt * rho*dx*dy * sum(vdirf(ebflag_macy));
            % % xdirf = 0.5 * (udirf(1:nx+2,:) + udirf(2:nx+3,:));
            % % ydirf = 0.5 * (vdirf(:,1:ny+2) + vdirf(:,2:ny+3));
            % % hydro_fx = -1/dt * rho*dx*dy * sum(xdirf(ebflag));
            % % hydro_fy = -1/dt * rho*dx*dy * sum(ydirf(ebflag));
            % hydro_fxdf = -1/dt * rho*dx*dy * sum(udirf(ebflag_macx));
            % hydro_fydf = -1/dt * rho*dx*dy * sum(vdirf(ebflag_macy));
            
            
            % [diffu diffv] = VelocityDiffusion(uold,vold,nx,ny,dx,dy,dt);
            % [gpx,gpy] = VelocityCorrector(pold,rho,nx,ny,dx,dy,dt);
            % xdirf = rho * ebvof_macx .* (gpx + diffu);
            % ydirf = rho * ebvof_macy .* (gpy + diffv);
            % % hydro_fx1 = dx*dy * sum(xdirf(ebflag_macx));
            % % hydro_fy1 = dx*dy * sum(ydirf(ebflag_macy));
            % xdirf = 0.5 * (xdirf(1:nx+2,:)+xdirf(2:nx+3,:));
            % ydirf = 0.5 * (ydirf(:,1:ny+2)+ydirf(:,2:ny+3));
            % hydro_fx = dx*dy * sum(xdirf(ebflag));
            % hydro_fy = dx*dy * sum(ydirf(ebflag));
            % % hydro_fx2 = dx*dy * sum(xdirf(:));
            
            % % hydro_force(end+1,1:3) = [time, hydro_fx, hydro_fy];
            % % plot(hydro_force(:,1),hydro_force(:,2), hydro_force(:,1),hydro_force(:,3));
            % % legend('Drag (Fx)', 'Lift (Fy)');
            % % hydro_force(end+1,1:5) = [time,hydro_fx,hydro_fy,hydro_fx1,hydro_fx2];
            % % plot(hydro_force(:,1),hydro_force(:,2), hydro_force(:,1),hydro_force(:,3), ...
                % % hydro_force(:,1),hydro_force(:,4), hydro_force(:,1),hydro_force(:,5));
            % % legend('Fx','Fy','Fx-MAC','Fx-total');
            % % xlim([0 max_time]);
            % hydro_force(end+1,1:5) = [time, hydro_fx, hydro_fy, hydro_fxdf,hydro_fydf];
            % plot(hydro_force(:,1),hydro_force(:,2), hydro_force(:,1),hydro_force(:,3), ...
                % hydro_force(:,1),hydro_force(:,4));
            % legend('Drag (Fx)', 'Lift (Fy)', 'Drag (direct)');
            % title(['hydrodynamics force ', num2str(hydro_fx)]);
            
            % subplot(subnrow,subncol,6);
            % plot(hydro_force(:,1),hydro_force(:,2)/(0.5*rho*UIn^2*cylinder_D), ...
                % hydro_force(:,1),hydro_force(:,4)/(0.5*rho*UIn^2*cylinder_D), ...
                % hydro_force(:,1),hydro_force(:,3)/(0.5*rho*UIn^2*cylinder_D));
            % legend('Cd (stress)', 'Cd (direct)', 'Cl (stress)');
            % title(['Cd=', num2str(hydro_fx/(0.5*rho*UIn^2*cylinder_D)), ...
            % ';Cl=', num2str(hydro_fy/(0.5*rho*UIn^2*cylinder_D))]);
            
            % update head-up plots
            drawnow;
        end
    end
end


if (1)
    jaxis = ny/2+2;
    
    figure;
    subplot(1,2,1);
    plot(validxs,vmac(2:nx+1,jaxis), validxs,vana(2:nx+1,jaxis));
    legend('sim','ana');
    title('vel');
    
    subplot(1,2,2);
    plot(validxs,0.5*(pres(2:nx+1,jaxis-1)+pres(2:nx+1,jaxis)), ...
    validxs,0.5*(pp(2:nx+1,jaxis-1)+pp(2:nx+1,jaxis)));
    legend('sim','ana');
    title('pres');
end




