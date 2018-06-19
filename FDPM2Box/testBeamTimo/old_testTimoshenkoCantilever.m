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

%% The old code


clear;

rho0 = 1; % we don't really need the density

E0 = 21.1 * 1e6;
nu0 = 0.3;

P0 = 10 * 1e3;
L0 = 10;
c0 = 1;

% gravity = [0; -1];
gravity = [0; 0];


% Timoshenko's cantilever solution
% we assume plain-strain state (which is different from plain-stress commonly used)
timo = refsol_TimoshenkoCantilever(E0,nu0, L0,c0,P0);


% setup neo-Hookean material
matl = materialNeoHookeanPlaneStrain(E0,nu0);


% particle size 
refine = 4;
% refine = 8;
nx = refine*10 + 5;
ny = refine*2 + 1;
h0 = L0 / nx;

% influence radius, if use p(2), must > 2
dilation = 2.1;
re = h0 * dilation;

% finite-increment-gradient stabilization
% fig_stab_alpha = 0.0;
fig_stab_alpha = 0.5;


% point number
numNodes = nx * ny;
% dof is (udisp,vdisp)
numDofs = numNodes * 2;

% generate points
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

% calculate Timoshenko's analytical solution
timosol_u = timo.ux(nodeX,nodeY);
timosol_v = timo.uy(nodeX,nodeY);
timosol_sxx = timo.sxx(nodeX,nodeY);
timosol_syy = timo.syy(nodeX,nodeY);
timosol_sxy = timo.sxy(nodeX,nodeY);
timosol_x = nodeX + timosol_u;
timosol_y = nodeY + timosol_v;

% this is tricky
% we would like to compare the converge, so the exact solution is used as BC
% dispBCDofs = reshape([fixed_nodes*2-1; fixed_nodes*2], [],1);
dispBCDofs = fdpmNode2Dof(fixed_nodes, 'xy');
dispBCVals = reshape([timosol_u(fixed_nodes), timosol_v(fixed_nodes)]', [],1);

% tracBCDofs = reshape([loaded_nodes*2-1; loaded_nodes*2], [],1);
tracBCDofs = fdpmNode2Dof(loaded_nodes, 'xy');
tracBCVals = reshape([zeros(length(loaded_nodes),1), ...
	nodeSize(loaded_nodes,2).*timosol_sxy(loaded_nodes)]', [],1);

% 
figure;
fdpmPlotNodes;

% 
disp('Begin build connnectivity');
tic;
clear conn;
% dNX = sparse(numNodes,numNodes);
% dNY = sparse(numNodes,numNodes);
% dNXX = sparse(numNodes,numNodes);
% dNXY = sparse(numNodes,numNodes);
% dNYY = sparse(numNodes,numNodes);
for i = 1:numNodes
	
	[neigh,rs,rx,ry,cutoff] = fdpmNeighborhood(i,nodePos,re);
	neigh2 = [neigh; i]; % particle included in the neighbourhood
	
	% weight function
	w = fdpmWeightFunc(rs(neigh),cutoff);
	
	% derivative shape function
	dN = fdpmShapeDer(rx(neigh),ry(neigh),w);
	
	% finite-increment-gradient correction to 1st-derivatives
	if (fig_stab_alpha > 0)
		hX = fig_stab_alpha*nodeSize(i,1);
		hY = fig_stab_alpha*nodeSize(i,2);
		
		[dN(1,:),dN(2,:)] = fdpmShapeDerFIGStab(dN(1,:),dN(2,:),dN(3,:),dN(4,:),dN(5,:), hX,hY);
	end
	
	% save connection
	% conn(i,neigh) = 1;
	% conn(i,i) = 1;
	conn(i).neigh = neigh;
	conn(i).neigh2 = neigh2;
	conn(i).dNX = dN(1,:);
	conn(i).dNY = dN(2,:);
	conn(i).dNXX = dN(3,:);
	conn(i).dNXY = dN(4,:);
	conn(i).dNYY = dN(5,:);
	
	
	% dNX(i,neigh2) = dN(1,:);
	% dNY(i,neigh2) = dN(2,:);
	% dNXX(i,neigh2) = dN(3,:);
	% dNXY(i,neigh2) = dN(4,:);
	% dNYY(i,neigh2) = dN(5,:);
end
% transpose
for i = 1:numNodes
	ineigh = conn(i).neigh2(:)';
	nneigh = numel(ineigh);
	
	conn(i).dNXt = zeros(1,nneigh);
	conn(i).dNYt = zeros(1,nneigh);
	
	for jj = 1:nneigh
		j = ineigh(jj);
		
		%
		ii = find(conn(j).neigh2 == i);
		
		conn(i).dNXt(jj) = conn(j).dNX(ii);
		conn(i).dNYt(jj) = conn(j).dNY(ii);
	end
end
toc;
disp('End build connnectivity');


% allocate buffers
nodeF = zeros(2,2,numNodes);
nodeP = zeros(2,2,numNodes);
nodeS = zeros(2,2,numNodes);
nodeSigma = zeros(2,2,numNodes);





% incremental load step
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

% incremental load loop
for step = 1:maxStep
	%
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
		
		disp('Build stiffness matrix');
		tic;
		% build stiffness matrix
		K = sparse(numDofs,numDofs);
		% K = zeros(numDofs);
		
		for i = 1:numNodes
			% neighbor list
			ineigh = conn(i).neigh2(:)';
			nneigh = numel(ineigh);
			
			ineighpos = pos(ineigh,:);
			ineighdof = fdpmNode2Dof(ineigh, 'xy');
			
			% d/dX
			% DN = [dNX(i,:); dNY(i,:)];
			DN = [conn(i).dNX; conn(i).dNY];
						
			% deformation gradient
			Fi = ineighpos' * DN';
			Ji = det(Fi);
			
			%
			Voli = nodeVol(i);
			Voli = Voli * Ji;
			
			% Cauchy stress, neo-Hookean
			sigma = matl.sigma(Fi,Ji);
			
			% constitution matrix, neo-Hookean
			Dmat = matl.D(Fi,Ji);
			
			% d/dx
			% dN = sparse(Fi' \ DN);
			dN = Fi' \ DN;
			dNx = dN(1,:);
			dNy = dN(2,:);
			
			% pad0 = zeros(1,numNodes);
			pad0 = zeros(1,nneigh);
			
			% constitutive portion
			Bmat = [reshape([dNx; pad0], 1,[]); ...
				reshape([pad0; dNy], 1,[]); ...
				reshape([dNy; dNx], 1,[])];
			% Kc = nodeVol(i) * (Bmat' * sparse(Dmat) * Bmat);
			Kc = Voli * (Bmat' * Dmat * Bmat);
			
			% initial stress portion
			Gmat = [reshape([dNx; pad0], 1,[]); ...
				reshape([dNy; pad0], 1,[]); ...
				reshape([pad0; dNx], 1,[]); ...
				reshape([pad0; dNy], 1,[])];
			smat = [sigma, zeros(2,2); zeros(2,2),sigma];
			% Ks = nodeVol(i) * (Gmat' * sparse(smat) * Gmat);
			Ks = Voli * (Gmat' * smat * Gmat);
			
			% assemble
			% K = K + Kc + Ks;
			K(ineighdof,ineighdof) = K(ineighdof,ineighdof) + Kc + Ks;
		end
		
		K = sparse(K);
		
		rhs = -residual;
		
		% apply BC
		K(dispBCDofs,:) = 0;
		K(:,dispBCDofs) = 0;
		spI = speye(numDofs);
		K(dispBCDofs,:) = spI(dispBCDofs,:);
		rhs(dispBCDofs) = 0;
		toc;
		
		disp('Solve matrix');
		% matrix solver
		tic;
		sol = K \ rhs;
		relax = 1;
		pos = pos + relax * reshape(sol,2,numNodes)';
		toc;
		
		% find real nodal force
		for i = 1:numNodes
			% neighbor list
			ineigh = conn(i).neigh2(:)';
			% neighbor positions
			ineighpos = pos(ineigh,:);
			
			% d/dX
			DN = [conn(i).dNX; conn(i).dNY];
			
			% deformation gradient
			Fi = ineighpos' * DN';
			Ji = det(Fi);
			
			% neo-Hookean
			Si = matl.S(Fi,Ji);
			Pi = Fi * Si;
			sigma = matl.sigma(Fi,Ji);
			nodeS(:,:,i) = Si;
			nodeP(:,:,i) = Pi;
			nodeSigma(:,:,i) = sigma;
		end
		
		% 
		% for i = 1:numNodes
			% ineigh = conn(i).neigh2;
			% pvec = reshape(nodeP(:,:,ineigh), 2,[]);
			% volvec = nodeVol(ineigh)';
			% dervec = bsxfun(@times, [conn(i).dNXt; conn(i).dNYt], volvec);
			
			% Ti = pvec * reshape(dervec,[],1);
			% T(:,i) = Ti;
		% end
		% internal force T
		T = zeros(2,numNodes);
		for i = 1:numNodes
			ineigh = conn(i).neigh2;
			
			DN = [conn(i).dNX; conn(i).dNY];
			
			Pi = nodeP(:,:,i);
			
			% force by i to neighbors j
			Ti = Pi * DN;
			Ti = Ti .* nodeVol(i);
			
			T(:,ineigh) = T(:,ineigh) + Ti;
		end
		T = reshape(T, [],1);
		
		% 
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







