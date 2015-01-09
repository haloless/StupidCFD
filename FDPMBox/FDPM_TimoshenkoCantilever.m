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

% ## FDPM_cantilever

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-31

clc;
clear all;

rho0 = 1; % we don't really need the density

% E0 = 1e4;
% nu0 = 0.25;
% P0 = 1;
% L0 = 24;
% c0 = 2;
E0 = 21.1 * 1e6;
nu0 = 0.3;
P0 = 10 * 1e3;
L0 = 10;
c0 = 1;

d0 = c0 * 2;
I0 = d0^3 / 12;

% Lame's constants
lambda0 = E0*nu0 / (1+nu0) / (1-2*nu0);
mu0 = E0 / 2 / (1+nu0);

% Timoshenko's cantilever
% we assume plain-strain state (which is different from plain-stress commonly used)
Ebar = E0 / (1-nu0^2);
nubar = nu0 / (1-nu0);
timo_ux = @(x,y) -P0/(6*Ebar*I0).* y .*((6*L0-3*x).*x + (2+nubar).*(y.^2-c0^2));
timo_uy = @(x,y) P0/(6*Ebar*I0).*(3*nubar*y.^2.*(L0-x) + (4+5*nubar)*c0^2*x + (3*L0-x).*x.^2);
timo_sxx = @(x,y) P0/I0 * (L0-x).*y;
timo_syy = @(x,y) zeros(size(x));
timo_sxy = @(x,y) P0/(2*I0) * (c0^2 - y.^2);


refine = 4;
nx = refine*10 + 5;
ny = refine*2 + 1;
h0 = L0 / nx;
% if use p(2), must > 2
dilation = 2.1;
re = h0 * dilation;

weightFunc = @(q) (q<1) .* (1 - 6*q.^2 + 8*q.^3 - 3*q.^4);

% finite-increment-gradient stabilization
do_FIGStab = 0;
% Taylor-series-expansion based stabilization
do_TEBStab = 0;
if(length(find([do_FIGStab do_TEBStab])) > 1)
    error('Only one stabilization is allowed.');
end


numNodes = nx * ny;
[gridX,gridY] = ndgrid(linspace(h0/2,L0-h0/2,nx),linspace(-c0+h0/2,c0-h0/2,ny));
nodeX = reshape(gridX,[],1);
nodeY = reshape(gridY,[],1);
nodePos = [nodeX,nodeY];
nodeVol = h0^2 * ones(numNodes,1);
nodeMass = rho0 * nodeVol;
nodeSize = h0 * ones(numNodes,2);
% BC
fixed_nodes = [1:nx:numNodes];
loaded_nodes = [nx:nx:numNodes];

if (0)
    figure;
    plot(nodeX,nodeY,'o', ...
    nodeX(fixed_nodes),nodeY(fixed_nodes),'x', ...
    nodeX(loaded_nodes),nodeY(loaded_nodes),'s');
    legend('node','disp','trac');
    axis([0 L0 -c0 c0]);
    axis('equal');
    % return
end

conn = sparse(numNodes,numNodes);
dNX = sparse(numNodes,numNodes);
dNY = sparse(numNodes,numNodes);
dNXX = sparse(numNodes,numNodes);
dNXY = sparse(numNodes,numNodes);
dNYY = sparse(numNodes,numNodes);
for i = 1:numNodes
    % xc = nodePos(i,1); yc = nodePos(i,2);
    % rx = nodeX - xc;
    % ry = nodeY - yc;
    % rs = sqrt(rx.^2 + ry.^2);
    % rs(i) = Inf; % exclude i-th particle itself
    % neigh = find(rs < re); % particle itself not in the neighbourhood
    
    [neigh,rs,rx,ry,cutoff] = FDPM_Neighborhood(i,nodePos,re);
    neigh2 = [neigh; i]; % particle included in the neighbourhood
    
    completeness = 2;
    switch completeness
    case {1}
        P = [rx(neigh),ry(neigh)];
    case {2}
        P = [rx(neigh),ry(neigh), 1/2*rx(neigh).^2,rx(neigh).*ry(neigh),1/2*ry(neigh).^2];
    otherwise
        error('Unsupported completeness %d',completeness);
    end
    w = weightFunc(rs(neigh)./cutoff);
    vol = nodeVol(neigh);
    
    Amat = P' * diag(w.*vol) * P;
    dN = (Amat' \ P') * diag(w.*vol);
    
    % save connection
    conn(i,neigh) = 1;
    conn(i,i) = 1;
    % save shape function
    dNX(i,neigh) = dN(1,:);
    dNY(i,neigh) = dN(2,:);
    dNX(i,i) = -sum(dNX(i,neigh));
    dNY(i,i) = -sum(dNY(i,neigh));
    
    dNXX(i,neigh) = dN(3,:);
    dNXX(i,i) = -sum(dNXX(i,neigh));
    dNXY(i,neigh) = dN(4,:);
    dNXY(i,i) = -sum(dNXY(i,neigh));
    dNYY(i,neigh) = dN(5,:);
    dNYY(i,i) = -sum(dNYY(i,neigh));
    
    % finite-increment-gradient
    if (do_FIGStab)
        alpha = 0.5;
        hX = alpha*nodeSize(i,1); hY = alpha*nodeSize(i,2);
        % dNX(i,neigh) = dNX(i,neigh) + hX*dN(3,:) + hY*dN(4,:);
        % dNY(i,neigh) = dNY(i,neigh) + hY*dN(5,:) + hX*dN(4,:);
        dNX(i,:) = dNX(i,:) + hX*dNXX(i,:) + hY*dNXY(i,:);
        dNY(i,:) = dNY(i,:) + hY*dNYY(i,:) + hX*dNXY(i,:);
    end
end

if (0) % check
    for i = 1:numNodes
        tx = eye(2,2);
        u = nodePos * tx';
        Fx = u' * [dNX(i,:)',dNY(i,:)'];
        if (norm(Fx-tx) > 1e-8)
            disp(['check1: node ', int2str(i), ' dXdX error']);
            disp([Fx tx]);
        end
        
        % nodex = nodeX * 1 + nodeY * 2;
        % nodey = nodeX * 3 + nodeY * 4;
        tx = [1 2;3 4];
        u = nodePos * tx';
        Fx = u' * [dNX(i,:)',dNY(i,:)'];
        if (norm(Fx-tx) > 1e-8)
            disp(['check2: node ', int2str(i), ' dxdX error']);
            disp([Fx tx]);
        end
        
        % this is for order(2) only
        u = nodePos.^2;
        Fx = u' * [dNX(i,:)',dNY(i,:)'];
        tx = [2*nodePos(i,1),0; 0,2*nodePos(i,2)];
        if (norm(Fx-tx) > 1e-8)
            disp(['check3: node ', int2str(i), ' dxdX error']);
            disp([Fx tx]);
        end
        
        u = [1/2*nodePos(:,1).^2, nodePos(:,1).*nodePos(:,2), 1/2*nodePos(:,2).^2];
        Fx = u' * [dNXX(i,:)', dNXY(i,:)', dNYY(i,:)'];
        tx = [1 0 0; 0 1 0; 0 0 1];
        if (norm(Fx-tx) > 1e-8)
            disp(['check4: node ', int2str(i), ' dxdX error']);
            disp([Fx tx]);
        end
    end
end

position = nodePos;
velocity = zeros(numNodes,2);
% gravity = [0; -1];
gravity = [0; 0];

nodeF = zeros(2,2,numNodes);
nodeP = zeros(2,2,numNodes);
nodeS = zeros(2,2,numNodes);
nodeSigma = zeros(2,2,numNodes);

if (1) % static load
    timosol_u = timo_ux(nodeX,nodeY);
    timosol_v = timo_uy(nodeX,nodeY);
    timosol_sxx = timo_sxx(nodeX,nodeY);
    timosol_syy = timo_syy(nodeX,nodeY);
    timosol_sxy = timo_sxy(nodeX,nodeY);
    timosol_x = nodeX + timosol_u;
    timosol_y = nodeY + timosol_v;
    if (0)
        figure;
        plot(timosol_x,timosol_y,'o');
        axis([0 L0 -c0 c0]);
        axis equal; 
        % return
    end
    
    numDofs = numNodes * 2;
    
    dispBCDofs = reshape([fixed_nodes*2-1; fixed_nodes*2], [],1);
    dispBCVals = reshape([timosol_u(fixed_nodes), timosol_v(fixed_nodes)]', [],1);
    
    tracBCDofs = reshape([loaded_nodes*2-1; loaded_nodes*2], [],1);
    tracBCVals = reshape([zeros(length(loaded_nodes),1), ...
        nodeSize(loaded_nodes,2).*timosol_sxy(loaded_nodes)]', [],1);
    
    maxStep = 1;
    % totalF = reshape([nodeMass'*gravity(1);nodeMass'*gravity(2)],[],1);
    totalF = zeros(numDofs,1);
    totalF(tracBCDofs) = tracBCVals;
    incrF = totalF / maxStep;
    
    totalU = zeros(numDofs,1);
    totalU(dispBCDofs) = dispBCVals;
    incrU = totalU / maxStep;
    
    % pos = reshape(nodePos',[],1);
    pos = nodePos;
    loadF = zeros(numDofs,1);
    residual = zeros(numDofs,1);
    
    for step = 1:maxStep
        pos = pos + reshape(incrU,2,numNodes)';
        loadF = loadF + incrF;
        residual = residual - incrF;
        
        disp(['load increment ', int2str(step)]);
        
        tol_abs = -1;
        tol_rel = 1e-9;
        maxIter = 50;
        converged = 0;
        
        % Newton-Raphson inner loop
        for iter = 1:maxIter
            K = sparse(numDofs,numDofs);
            % build stiffness matrix
            for i = 1:numNodes
                % d/dX
                DN = [dNX(i,:); dNY(i,:)];
                % deformation gradient
                Fi = pos' * DN';
                Ji = det(Fi);
                % Cauchy stress, neo-Hookean
                sigma = mu0/Ji*(Fi*Fi'-eye(2)) + lambda0*log(Ji)/Ji*eye(2);
                % constitution matrix, neo-Hookean
                lambdai = lambda0 / Ji;
                mui = (mu0-lambda0*log(Ji)) / Ji;
                Dmat = [lambdai+2*mui, lambdai, 0; ...
                    lambdai, lambdai+2*mui, 0; ...
                    0, 0, mui];
                
                % d/dx
                dN = sparse(Fi' \ DN);
                dNx = dN(1,:);
                dNy = dN(2,:);
                
                pad0 = zeros(1,numNodes);
                % constitutive portion
                Bmat = [reshape([dNx; pad0], 1,[]); ...
                    reshape([pad0; dNy], 1,[]); ...
                    reshape([dNy; dNx], 1,[])];
                Kc = nodeVol(i) * (Bmat' * sparse(Dmat) * Bmat);
                
                % initial stress portion
                Gmat = [reshape([dNx; pad0], 1,[]); ...
                    reshape([dNy; pad0], 1,[]); ...
                    reshape([pad0; dNx], 1,[]); ...
                    reshape([pad0; dNy], 1,[])];
                smat = [sigma, zeros(2,2); zeros(2,2),sigma];
                Ks = nodeVol(i) * (Gmat' * sparse(smat) * Gmat);
                
                % assemble
                K = K + Kc + Ks;
            end            
            rhs = -residual;
            
            % apply BC
            K(dispBCDofs,:) = 0;
            K(:,dispBCDofs) = 0;
            spI = speye(numDofs);
            K(dispBCDofs,:) = spI(dispBCDofs,:);
            rhs(dispBCDofs) = 0;
            
            u = K \ rhs;
            relax = 1;
            pos = pos + relax * reshape(u,2,numNodes)';
            
            % find real nodal force
            T = zeros(2,numNodes);
            nodeP = zeros(2,2,numNodes);
            for i = 1:numNodes
                Fi = pos' * [dNX(i,:); dNY(i,:)]';
                % neo-Hookean
                J = det(Fi);
                invCi = inv(Fi' * Fi);
                Si = mu0*(eye(2)-invCi) + lambda0*log(J)*invCi;
                Pi = Fi * Si;
                nodeS(:,:,i) = Si;
                nodeP(:,:,i) = Pi;
                
                sigma = mu0/J*(Fi*Fi'-eye(2)) + lambda0*log(J)/J*eye(2);
                nodeSigma(:,:,i) = sigma;
            end
            nodeP = reshape(nodeP,2,[]);
            nodeV = diag(reshape([nodeVol'; nodeVol'], [],1));
            for i = 1:numNodes
                Ti = nodeP * nodeV * reshape([dNX(:,i)'; dNY(:,i)'], [],1);
                T(:,i) = Ti;
            end
            T = reshape(T, [],1);
            
            residual = T - loadF;
            % residual at fixed nodes is nonsense
            residual(dispBCDofs) = 0;
            
            if (1)
                plot(pos(:,1),pos(:,2),'o');
                axis equal;
                title(['outer:',int2str(step),'/inner:',int2str(iter)]);
                drawnow;
            end
            
            normRes = norm(residual);
            normLoad = norm(loadF);
            disp(['Newton iteration=',int2str(iter),'; |res|/|load|=',num2str(normRes/(normLoad+eps))]);
            if (normRes<=tol_rel*normLoad || normRes<=tol_abs)
                converged = 1;
                break;
            end
        end
        
        if (converged)
            disp(['load step ',int2str(step),' converged']);
        else
            warning('load step %d did NOT converged to wanted threshold', step);
        end
    end % load step end
    
    % post-processing
    
    figure;
    xpos = reshape(pos(:,1), nx,ny);
    ypos = reshape(pos(:,2), nx,ny);
    xres = reshape(residual, 2,nx*ny); xres = xres(1,:); xres = reshape(xres, nx,ny);
    yres = reshape(residual, 2,nx*ny); yres = yres(2,:); yres = reshape(yres, nx,ny);
    surf(xpos',ypos',sqrt(xres.^2+yres.^2)');
    title('residual of balance equation');
    
    figure;
    plot(xpos(:),ypos(:),'o', nodeX+timosol_u,nodeY+timosol_v, 'x');
    legend('FDPM','Timo');
    axis([0 L0 -c0 +c0]); axis equal;
    
    figure;
    PK2st = reshape(nodeSigma, 4,nx*ny);
    PK2st(:,fixed_nodes) = NaN;
    S11 = reshape(PK2st(1,:), nx,ny);
    S21 = reshape(PK2st(2,:), nx,ny);
    S12 = reshape(PK2st(3,:), nx,ny);
    S22 = reshape(PK2st(4,:), nx,ny);
    subplot(2,2,1);
    surf(xpos',ypos', S11');
    title('S11'); colorbar; view(0,90);
    subplot(2,2,2);
    surf(xpos',ypos', S21');
    title('S21'); colorbar; view(0,90);
    subplot(2,2,3);
    surf(xpos',ypos', S12');
    title('S12'); colorbar; view(0,90);
    subplot(2,2,4);
    surf(xpos',ypos', S22');
    title('S22'); colorbar; view(0,90);
    
    figure;
    proby = sub2ind([nx,ny], 1:nx, repmat((ny+1)/2,1,nx));
    plot(nodeX(proby),S11(proby), nodeX(proby),S21(proby), ...
        nodeX(proby),S12(proby), nodeX(proby),S22(proby), ...
        nodeX(proby),timosol_sxx(proby), ...
        nodeX(proby),timosol_sxy(proby), ...
        nodeX(proby),timosol_syy(proby));
    legend('S11','S21','S12','S22', 'sxx-timo', 'sxy-timo', 'syy-timo');
    
    figure;
    plot(xpos(proby),ypos(proby), timosol_x(proby),timosol_y(proby));
    legend('FDPM', 'Timo');
    
    if (1) % do some reporting
        disp(['num=',int2str(numNodes), ...
            '; nx=', int2str(nx), '; ny=', int2str(ny), ...
            '; dh=', num2str(h0), '; re=', num2str(dilation)]);
        
        errL2 = sum((xpos(:)-nodeX-timosol_u).^2) + sum((ypos(:)-nodeY-timosol_v).^2);
        errL2 = sqrt(errL2 * h0^2);
        disp(['n=', int2str(numNodes), '; dh=', num2str(h0), '; L2=', num2str(errL2)]);
    end

end % static load






