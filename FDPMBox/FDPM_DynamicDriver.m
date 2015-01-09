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

% ## FDPM_DynamicDriver

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-08-20

% function [ ret ] = FDPM_DynamicDriver ()

if (1) % dynamic run
    dt = 5/2 * 1e-3;
    dt = 1 * 1e-3;
    dt = 2.5 * 1e-4;
    time = 0;
    energy = [];
    times = [];
    
    for step = 1:5000
        time = time + dt;
        
        nodeJ = zeros(numNodes,1);
        for i = 1:numNodes
            Fi = position' * [dNX(i,:)',dNY(i,:)'];
            
            % % St. Venant-Kirchhoff
            % Ei = 1/2 * (Fi'*Fi - eye(2));
            % Si = 2*mu0*Ei + lambda0*trace(Ei)*eye(2);
            
            % neo-Hookean
            J = det(Fi);
            Ci = Fi' * Fi;
            invCi = inv(Ci);
            if(0)
                Si = mu0*(eye(2)-invCi) + lambda0*log(J)*invCi;
            else
                % bi = Fi * Fi';
                % sig = mu0 * J^(-5/3) * (bi - trace(bi)/3*eye(2));
                % pres = lambda0*(J-1);
                % sig = sig + pres*eye(2);
                % Si = J * inv(Fi) * (sig + pres*eye(2)) * inv(Fi');
                Ib = trace(Ci);
                pres = lambda0 * (J-1);
                Si = mu0*J^(-2/3)*(eye(2)-Ib/3*invCi) + J*pres*invCi;
            end
            
            
            Pi = Fi * Si;
            nodeF(:,:,i) = Fi;
            nodeP(:,:,i) = Pi;
            
            nodeJ(i) = J;
        end
        
        Ps = reshape(nodeP,2,[]);
        Vs = diag(reshape([nodeVol';nodeVol'],[],1));
        
        xnew = position;
        vnew = velocity;
        acc = zeros(numNodes,2);
        
        for i = 1:numNodes
            Ti = Ps * Vs * reshape([dNX(:,i)';dNY(:,i)'],[],1);
            ai = -Ti/nodeMass(i) + gravity;
            acc(i,:) = ai';
        end
        
        % if step==1
            % vnew = vnew + 1/4*dt * acc;
            % xnew = xnew + 1/2*dt * vnew;
            % vnew = vnew + 3/4*dt * acc;
        % else
            % vnew = vnew + dt*acc;
            % xnew = xnew + dt*vnew;
        % end
        vnew = vnew + dt*acc;
        xnew = xnew + dt*vnew;
        
        % fixed BC
        vnew(fixed_nodes,:) = 0;
        xnew(fixed_nodes,:) = position(fixed_nodes,:);
        
        velocity = vnew;
        position = xnew;
        
        
        if (mod(step,10) == 0)
            prompt = ['step=',int2str(step),';time=',num2str(time), ...
            ';Jmax=',num2str(max(nodeJ)), ';Jmin=',num2str(min(nodeJ))];
            % plot(position(:,1),position(:,2),'o');
            plot(position(:,1),position(:,2),'o', nodeX,nodeY,'x');
            title(prompt);
            axis('equal');
            axis([0 L0*1.2 -c0*4 c0*4]);
            hold on;
            
            quiver(position(:,1),position(:,2),velocity(:,1),velocity(:,2));
            hold off;
            
            drawnow;
            
            times(end+1) = time;
            energy(end+1) = 1/2 * sum(nodeMass .* (velocity(:,1).^2+velocity(:,2).^2));
        end
    end
    if (1)
        figure;
        plot(times,energy);
        title('Ek');
    end
end % dynamic run

% return
% end
