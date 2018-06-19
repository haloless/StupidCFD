
%% old code


clear;




E0 = 21.1 * 1e6;
nu0 = 0.3;

P0 = 10 * 1e3;
L0 = 10;
c0 = 1;

scale = 1.0e4;
E0 = E0 / scale;
P0 = P0 / scale;

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
nodeSize = h0 * ones(numNodes,2);

%
% BC
%
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
% displacement at fixed end
dispBCDofs = fdpmNodeDof(fixed_nodes, 'xy');
dispBCVals = utilZipVec(timosol_u(fixed_nodes), timosol_v(fixed_nodes));
% traction at free end
tracBCDofs = fdpmNodeDof(loaded_nodes, 'xy');
tracBCVals = utilZipVec([], nodeSize(loaded_nodes,2).*timosol_sxy(loaded_nodes));

if (1)
    % plot fdpm node mesh
    figure;
    fdpmPlotNodes;
end


%% build connnectivity

disp('Begin build connnectivity');
tic;
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
    conn(i).numNeigh = numel(neigh2);
	conn(i).neigh2 = neigh2'; % neighbor list is a row vector
	conn(i).dNX = dN(1,:);
	conn(i).dNY = dN(2,:);
	conn(i).dNXX = dN(3,:);
	conn(i).dNXY = dN(4,:);
	conn(i).dNYY = dN(5,:);
end
toc;
disp('End build connnectivity');


%% allocate buffers

% [fx,fy]
fint = zeros(numDofs,1);
fresid = zeros(numDofs,1);

% displacement [u,v]
uv = zeros(numDofs,1);

% elastic log strain [11,22,12]
nodeEpsE = zeros(3,numNodes);

% cauchy stress [11,22,12,33] 
nodeSigma = zeros(4,numNodes);

% deform grad
nodeF = zeros(3,3,numNodes);

% coord [x1,y1,x2,y2,...]
nodeCoord = nodePos';

% initialize
for i = 1:numNodes
    nodeF(:,:,i) = eye(3);
end

% save old values
uv_old = uv;
nodeEpsE_old = nodeEpsE;
nodeF_old = nodeF;

% create initial stiffness matrix
Ktan0 = sparse(numDofs,numDofs);
D = E0/(1+nu0)/(1-2*nu0) * [1-nu0,nu0,0;nu0,1-nu0,0;0,0,0.5-nu0];
for i = 1:numNodes
    ineigh = conn(i).neigh2;
    ineighdof = fdpmNodeDof(ineigh,'xy');
    ivol = nodeVol(i);
    
    B = fdpmFormBmat(conn(i).numNeigh, conn(i).dNX, conn(i).dNY);
    Ke = B' * D * B * ivol;
    
    Ktan0(ineighdof,ineighdof) = Ktan0(ineighdof,ineighdof) + Ke;
end

%% solve infinitesimal problem
if 1
    fext = zeros(numDofs,1);
    fext(tracBCDofs) = tracBCVals;
    
    frhs = fext - fint;
    
    sol = fdpmSolve(Ktan0,frhs, dispBCDofs,dispBCVals);
    
    uv = reshape(sol, 2,numNodes);
    
    nodeCoord = nodeCoord + uv;
    xpos = reshape(nodeCoord(1,:), nx,ny);
    ypos = reshape(nodeCoord(2,:), nx,ny);
    
    figure;
    proby = sub2ind([nx,ny], 1:nx, repmat((ny+1)/2,1,nx));
    plot(xpos(proby),ypos(proby), timosol_x(proby),timosol_y(proby));
    legend('FDPM', 'Timo');
end

% nodeF = zeros(2,2,numNodes);
% nodeP = zeros(2,2,numNodes);
% nodeS = zeros(2,2,numNodes);
% nodeSigma = zeros(2,2,numNodes);


return


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

if 0


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
			ineighdof = fdpmNodeDof(ineigh, 'xy');
			
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
		nodeP = zeros(2,2,numNodes);
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
		
		% % 
		% for i = 1:numNodes
			% ineigh = conn(i).neigh2;
			% pvec = reshape(nodeP(:,:,ineigh), 2,[]);
			% volvec = nodeVol(ineigh)';
			% dervec = bsxfun(@times, [conn(i).dNXt; conn(i).dNYt], volvec);
			
			% Ti = pvec * reshape(dervec,[],1);
			% T(:,i) = Ti;
		% end
		% T = reshape(T, [],1);
		
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

end

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







