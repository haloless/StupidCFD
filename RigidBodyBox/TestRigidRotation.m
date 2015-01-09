
clc;
clear all;


% generate particle system
npart = 0;

nrad = 4;
nphi = 12;
for i = 1:nrad
for j = 1:nphi
    rad = 1.0/nrad * i;
    phi = 2*pi/nphi * j;
    
    npart = npart + 1;
    xpart(npart) = rad * cos(phi);
    ypart(npart) = rad * sin(phi);
    zpart(npart) = 0.0;
end
end

if (0)
    scatter3(xpart,ypart,zpart);
end

% center position
xcen = sum(xpart) / npart;
ycen = sum(ypart) / npart;
zcen = sum(zpart) / npart;


% mass
mpart = ones(size(xcen));
% assume point, zero MOI
Jpart = zeros(3,3,npart);

xpart0 = xpart;
ypart0 = ypart;
zpart0 = zpart;



fig_main = figure;
set(gcf, 'name','RigidRotation');
btn_restart = uicontrol('position',[20 20 60 20],'style','toggle','string','restart');
btn_exit = uicontrol('position',[100 20 60 20],'style','toggle','string','exit');
edt_omegax = uicontrol('position',[180 20 60 20],'style','edit','string','0.0');
edt_omegay = uicontrol('position',[260 20 60 20],'style','edit','string','0.0');
edt_omegaz = uicontrol('position',[340 20 60 20],'style','edit','string','1.0');

fig_main_pos = get(fig_main,'Position');
set(fig_main,'Position',[fig_main_pos(1),fig_main_pos(2)-fig_main_pos(4),fig_main_pos(3),fig_main_pos(4)*2]);

% fig_err = figure;


% if (1)
while get(btn_exit,'value') == 0 % main loop
    set(btn_restart,'value',0);
    
    xpart = xpart0;
    ypart = ypart0;
    zpart = zpart0;
    J = DeriveShiftedMOI(xcen,ycen,zcen, npart, xpart,ypart,zpart, mpart,Jpart);
    omega = [ 0.0; 0.0; 1.0 ];
    % omega = [ 0.25; 0.5; 1.0 ];
    if (1)
        omega(1) = str2double(get(edt_omegax,'string'));
        omega(2) = str2double(get(edt_omegay,'string'));
        omega(3) = str2double(get(edt_omegaz,'string'));
    end
    
    L = J * omega;
    omega0 = omega;
    E = 0.5 * omega0' * L;
    L0 = L;
    E0 = E;
    
    e3 = L0 / norm(L0);
    e1 = [1; 0; 0];
    e1 = e1 - dot(e1,e3)*e3;
    e1 = e1 / norm(e1);
    e2 = [0; 1; 0];
    e2 = e2 - dot(e2,e3)*e3;
    e2 = e2 - dot(e2,e1)*e1;
    e2 = e2 / norm(e2);
    Pmat = [e1, e2, e3];
    
    time = 0;
    dt = 1.0e-3;
    max_step = 99999;
    
    Lplot = [time, L0', E0];
    
    for step = 1:max_step
        if get(btn_exit,'value')~=0, break, end
        if get(btn_restart,'value')~=0, break, end
        
        time = time + dt;
        
        if (0)
            gyro_mom = cross(omega, L);
            acc = J \ (-gyro_mom);
            
            % omega_adv = omega + 0.5*dt*acc;
            omega_adv = omega + dt*acc;
            omega = omega + dt*acc;
            
            [xpart,ypart,zpart] = UpdatePartPosition(xcen,ycen,zcen, ...
            omega_adv,dt,npart,xpart,ypart,zpart);
            
            %
            J = DeriveShiftedMOI(xcen,ycen,zcen, npart, xpart,ypart,zpart, mpart,Jpart);
            L = J * omega;
            E = 0.5 * omega' * L;
        end % Euler stepping
        
        if (1) % predictor-corrector
            %
            [xpart1,ypart1,zpart1] = UpdatePartPosition(xcen,ycen,zcen, ...
            omega,dt, npart,xpart,ypart,zpart);
            %
            J1 = DeriveShiftedMOI(xcen,ycen,zcen, npart, xpart1,ypart1,zpart1, mpart,Jpart);
            omega1 = J1 \ L;
            
            %
            omega_adv = 0.5 * (omega+omega1);
            [xpart,ypart,zpart] = UpdatePartPosition(xcen,ycen,zcen, ...
            omega_adv,dt,npart,xpart,ypart,zpart);
            %
            J = DeriveShiftedMOI(xcen,ycen,zcen, npart, xpart,ypart,zpart, mpart,Jpart);
            omega = J \ L;
            
            % this is not necessory
            L = J * omega;
            % energy
            E = 0.5 * omega' * L;
        end
        
        
        if mod(step,20) == 0 
            prompt = sprintf(['step=%d;time=%g;\n', ...
            'L0=(%g,%g,%g);L=(%g,%g,%g);\n', ...
            'w0=(%g,%g,%g);w=(%g,%g,%g);\n', ...
            'E0=%g;E=%g;'], ...
            step, time, ...
            L0(1),L0(2),L0(3), L(1),L(2),L(3), ...
            omega0(1),omega0(2),omega0(3), omega(1),omega(2),omega(3), ...
            E0,E);
            disp(prompt);
            
            % save for error plot
            Lplot(end+1,:) = [time, L', 0.5*omega'*L];
            
            % 
            figure(fig_main);
            plot_tx = 0;
            %plot_tx = 1;
            
            subplot(8,1,1:4);
            if plot_tx
                xyz_tx = inv(Pmat) * [xpart; ypart; zpart];
                scatter3(xyz_tx(1,:), xyz_tx(2,:), xyz_tx(3,:));
            else
                scatter3(xpart,ypart,zpart);
            end
            hold on;
            
            scatter3(xcen,ycen,zcen,'*');
            
            vel = cross(repmat(omega,[1,npart]), [xpart-xcen;ypart-ycen;zpart-zcen], 1);
            if plot_tx
                vel = inv(Pmat) * vel;
                quiver3(xyz_tx(1,:), xyz_tx(2,:), xyz_tx(3,:), vel(1,:),vel(2,:),vel(3,:));
            else
                quiver3(xpart,ypart,zpart, vel(1,:),vel(2,:),vel(3,:));
            end
            hold off;
            
            % legend('Particle','Center','Velocity');
            title(prompt);
            axis equal;
            axis([-1,1, -1,1, -1,1]);
            
            subplot(8,1,5);
            plot(Lplot(:,1), Lplot(:,2)-L0(1));
            subplot(8,1,6);
            plot(Lplot(:,1), Lplot(:,3)-L0(2));
            subplot(8,1,7);
            plot(Lplot(:,1), Lplot(:,4)-L0(3));
            
            subplot(8,1,8);
            plot(Lplot(:,1), Lplot(:,5)-E0);
            
            drawnow;
            
        end
    end % end simulation loop
    
    % if (0) % plot angular momentum error
        % figure(fig_err);
        
        % subplot(3,1,1);
        % plot(Lplot(:,1), Lplot(:,2)-L0(1));
        % subplot(3,1,2);
        % plot(Lplot(:,1), Lplot(:,3)-L0(2));
        % subplot(3,1,3);
        % plot(Lplot(:,1), Lplot(:,4)-L0(3));
    % end
    
end % end main loop

if (0)
    close(fig_main);
    close(fig_err);
end


