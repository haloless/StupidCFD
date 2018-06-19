%%Timoshenko Cantilever
% Infinitesimal deformation
% Elastic solution

clear;


%% physical properties

% coordinate system
prob_type = 1; % plane-strain
% prob_type = 2; % axisymmetric
disp(['prob_type=',int2str(prob_type)]);

% Young, Poisson
E0 = 21.1 * 1e6;
% nu0 = 0.3;
nu0 = 0.49999;

P0 = 10 * 1e3; % end loading

L0 = 10; % length
c0 = 1; % half width

stiff_scale = 1.0e5;
E0 = E0 / stiff_scale;
P0 = P0 / stiff_scale;

% Timoshenko's cantilever solution
timo = refsol_TimoshenkoCantilever(E0,nu0, L0,c0,P0);

% D is the elastic tangent modulus
D = materialVonMises(zeros(4,1), E0,nu0, -1, prob_type); 
if prob_type == 1
    Dmat = D(1:3,1:3);
elseif prob_type == 2
    Dmat = D;
end

%% geometry

% nrefine = 3;
nrefine = 5;
nx = L0 * nrefine + 1;
ny = c0*2 * nrefine + 1;

% influence radius, if use p(2), must > 2
h0 = L0 / (nx-1);
% dilation = 1.2;
% dilation = 1.5;
dilation = 1.8;
% dilation = 2.1;
re = h0 * dilation;


%% generate points and boundary conditions

% generation routine should setup the following data
numNodes = 0;
% particle states
nodeX = [];
nodeY = [];
nodeVol = [];
% BC Dirichlet
dispBCDofs = [];
dispBCVals = [];
% BC loading
tracBCDofs = [];
tracBCVals = [];

disp('Begin generate particles & BC');
if prob_type == 1
    % plane-strain
    
    % 
    loaded_nodes = [];
    loaded_force = [];
    % 
    fixed_nodes_x = [];
    fixed_disp_x = [];
    fixed_nodes_y = [];
    fixed_disp_y = [];
    
    for j = 1:ny
    for i = 1:nx
        xx = 0 + (i-1)*h0;
        yy = 0 + (j-1)*h0 - c0;
        
        numNodes = numNodes + 1;
        nodeX(numNodes,1) = xx;
        nodeY(numNodes,1) = yy;
        
        lx = h0; ly = h0;
        if i == 1 || i == nx
            lx = lx * 0.5;
        end
        if j == 1 || j == ny
            ly = ly * 0.5;
        end
        nodeVol(numNodes,1) = lx * ly;
        
        if i == 1 % displacement BC
        % if i==1 && (j==1 || j==ny)
            fixed_nodes_y(end+1,1) = numNodes;
            fixed_nodes_x(end+1,1) = numNodes;
            fixed_disp_y(end+1,1) = timo.uy(xx,yy);
            fixed_disp_x(end+1,1) = timo.ux(xx,yy);
        end
        if i == nx % traction BC
            ss = timo.sxy(xx,yy) * ly;
            loaded_nodes(end+1,1) = numNodes;
            loaded_force(1:2,end+1) = [0, ss];
        end
        % if i==1
            % loaded_nodes(end+1,1) = numNodes;
            % loaded_force(1:2,end+1) = [-timo.sxx(xx,yy), -timo.sxy(xx,yy)] * ly;
        % end
    end
    end
    
    % setup Dirichlet BC
    tmp = fdpmNodeDof(fixed_nodes_y,'y');
    dispBCDofs = [dispBCDofs; tmp];
    dispBCVals = [dispBCVals; fixed_disp_y];
    tmp = fdpmNodeDof(fixed_nodes_x,'x');
    dispBCDofs = [dispBCDofs; tmp];
    dispBCVals = [dispBCVals; fixed_disp_x];
    
    % setup Loading
    tracBCDofs = fdpmNodeDof(loaded_nodes,'xy');
    tracBCVals = reshape(loaded_force, [],1);
    
else
    error('Unknown prob_type');
end
disp('End generate particles & BC');

nodePos = [nodeX,nodeY];

% dof is (udisp,vdisp)
numDofs = numNodes * 2;

if (1) % plot fdpm node mesh
    fdpmDriverPlotNodes;
end

if 1
    fdpmDriverBuildTri;
    nodeVol = nodeArea;
    % return;
end


%% build connnectivity

disp('Begin build connnectivity'); tic;
% fdpmDriverBuildConn;
% fdpmDriverBuildConnIntegBndry;
fdpmDriverBuildConnIntegBndry2;
% fdpmDriverBuildConnIntegDomain;
disp('End build connnectivity'); toc;
% return;

fdpmDriverBuildMoment;


%% allocate buffers

% [fx,fy]
fext = zeros(numDofs,1);
fint = zeros(numDofs,1);
frct = zeros(numDofs,1);

fstab = zeros(numDofs,1);

% displacement [u,v]
uv = zeros(numDofs,1);

% elastic log strain [11,22,12,33]
nodeEpsE = zeros(4,numNodes);

% cauchy stress [11,22,12,33] 
nodeSigma = zeros(4,numNodes);

% deform grad
nodeF = zeros(3,3,numNodes);

% coord [x1,y1,x2,y2,...]
nodeCoord = nodePos';
nodeCurr = nodeCoord;

% plastic
nodeAlpha = zeros(numNodes,1);
nodeEpbar = zeros(numNodes,1);

% nodeFlag = zeros(numNodes,1);

% initialize
for i = 1:numNodes
    nodeF(:,:,i) = eye(3);
end

% save old values
uv_old = uv;
nodeEpsE_old = nodeEpsE;
nodeF_old = nodeF;
nodeEpbar_old = nodeEpbar;
nodeAlpha_old = nodeAlpha;

%% create initial stiffness matrix
disp('Begin create initial stiffness matrix'); tic;

% assembler
Kassem = SpMatAssem(numDofs);


% stiffness
for i = 1:numNodes
    nneigh = conn(i).numNeigh;
    ineigh = conn(i).neigh2;
    ineighdof = fdpmNodeDof(ineigh,'xy');
    ivol = nodeVol(i);
    ineighx = nodeX(ineigh);
    ineighy = nodeY(ineigh);
    
    B = fdpmFormBmat(nneigh, conn(i).dNX,conn(i).dNY, prob_type, ineighx,ineighy);
    
    if prob_type == 1
        Ke = B' * D(1:3,1:3) * B * ivol;
    elseif prob_type == 2
        Ke = B' * D * B * ivol;
    end
    
    % assemble
    SpMatAssemBlockWithDof(Kassem, Ke,ineighdof,ineighdof);
end

if 0
    nodeFlag2 = nodeFlag;
    for i = 1:numTri
        if any(nodeFlag(tri(i,:)))
            % nodeFlag2(tri(i,:)) = 1;
        end
    end
    for i = 1:numNodes
        nneigh = conn(i).numNeigh;
        ineigh = conn(i).neigh2;
        ineighdof = fdpmNodeDof(ineigh,'xy');
        ivol = nodeVol(i);
        ineighx = nodeX(ineigh);
        ineighy = nodeY(ineigh);
        
        % if nodeFlag2(i), continue, end
        
        % TODO
        B = fdpmFormBmat(nneigh, conn(i).dNX,conn(i).dNY, prob_type, ineighx,ineighy);
        Bx = fdpmFormBmat(nneigh, conn(i).dNXX, conn(i).dNXY, prob_type, ineighx,ineighy);
        By = fdpmFormBmat(nneigh, conn(i).dNXY, conn(i).dNYY, prob_type, ineighx,ineighy);
        
        Kstab = Bx'*Dmat*Bx.*nodeMomXX(i) + By'*Dmat*By.*nodeMomYY(i);
        Kstab = Kstab + (Bx'*Dmat*By + By'*Dmat*Bx).*nodeMomXY(i);
        % Kstab = Kstab + (B'*Dmat*Bx + Bx'*Dmat*B).*nodeMomX(i);
        % Kstab = Kstab + (B'*Dmat*By + By'*Dmat*B).*nodeMomY(i);
        
        SpMatAssemBlockWithDof(Kassem, Kstab,ineighdof,ineighdof);
    end
end

% create sparse K
Ktan0 = SpMatCreateAndClear(Kassem);

disp('End create initial stiffness matrix'); toc;

% disp('Begin create initial stab matrix'); tic;
% fdpmDriverCalcStab;
% Kstab0 = Kstab;
% Ktan0 = Ktan0 + Kstab0;
% disp('End create initial stab matrix'); toc;

if 1
    fext(tracBCDofs) = tracBCVals;
    
    
    % residual force
    fres = frct + fext - fint - fstab;
    
    [dduv,dreact] = fdpmSolve(Ktan0,fres, dispBCDofs,dispBCVals);
    
    % update
    uv = uv + dduv;
    frct = frct + dreact;
    
    % current position
    nodeDisp = reshape(uv, 2,numNodes);
    nodeCurr = nodeCoord + nodeDisp;
    
    for i = 1:numNodes
        nneigh = conn(i).numNeigh;
        ineigh = conn(i).neigh2;
        ineighdof = fdpmNodeDof(ineigh,'xy');
        
        % updated coords
        % ineighpos = nodeCoord(:,ineigh) + reshape(uv(ineighdof), 2,nneigh);
        % ineighpos = nodeCurr(:,ineigh);
        ineighpos = nodeCoord(:,ineigh);
        
        % d/dX
        DN = [conn(i).dNX; conn(i).dNY];
        
        dNx = DN(1,:);
        dNy = DN(2,:);
        
        % current volume
        ivol = nodeVol(i);
        
        
        % gradient operator
        Bmat = fdpmFormBmat(nneigh, dNx,dNy, prob_type,ineighpos(1,:),ineighpos(2,:));
        % Gmat = fdpmFormGmat(nneigh, dNx,dNy, prob_type,ineighpos(1,:),ineighpos(2,:));
        
        epsE = Bmat * uv(ineighdof);
        
        sigma = D * epsE;
        
        if prob_type == 1
            Ti = Bmat.' * sigma(1:3) * ivol;
        elseif prob_type == 2
            Ti = Bmat.' * sigma * ivol;
        end
        
        % assemble
        fint(ineighdof) = fint(ineighdof) + Ti;
        
        % save values
        % nodeF(:,:,i) = Fi;
        % nodeEpsE(:,i) = epsE;
        % nodeEpbar(i) = epbar;
        % nodeAlpha(i) = alpha;
        % nodeSigma(:,i) = sigma;
        nodeSigma(1:3,i) = sigma(1:3);
        
    end
    
end




%% postprocessing

plot_disp_scale = 5.0;

% calculate Timoshenko's analytical solution
timosol_u = timo.ux(nodeX,nodeY);
timosol_v = timo.uy(nodeX,nodeY);
timosol_sxx = timo.sxx(nodeX,nodeY);
timosol_syy = timo.syy(nodeX,nodeY);
timosol_sxy = timo.sxy(nodeX,nodeY);
timosol_x = nodeX + timosol_u;
timosol_y = nodeY + timosol_v;

%
nodeDisp = nodeDisp.';

pos = nodeCurr.';
xpos = pos(:,1);
ypos = pos(:,2);

figure;
triplot(tri, nodeX+nodeDisp(:,1).*plot_disp_scale, nodeY+nodeDisp(:,2).*plot_disp_scale); 
hold on;
plot(nodeX+timosol_u.*plot_disp_scale,nodeY+timosol_v.*plot_disp_scale, 'rx');
hold off;
legend('FDPM','Timo');
% axis([0 L0 -c0 +c0]); 
axis equal;

proby = sub2ind([nx,ny], 1:nx, repmat((ny+1)/2,1,nx));

if 1
PK2st = reshape(nodeSigma, 4,nx*ny);
% PK2st(:,fixed_nodes) = NaN;
% S11 = reshape(PK2st(1,:), nx,ny);
% S21 = reshape(PK2st(2,:), nx,ny);
% S12 = reshape(PK2st(3,:), nx,ny);
% S22 = reshape(PK2st(4,:), nx,ny);
S11 = PK2st(1,:);
S21 = PK2st(2,:);
S12 = PK2st(3,:);
S22 = PK2st(4,:);

figure;
plot(nodeX(proby),S11(proby), nodeX(proby),S21(proby), ...
	nodeX(proby),S12(proby), nodeX(proby),S22(proby), ...
	nodeX(proby),timosol_sxx(proby), '.', ...
	nodeX(proby),timosol_sxy(proby), '.', ...
	nodeX(proby),timosol_syy(proby), '.');
legend('S11','S21','S12','S22', 'sxx-timo', 'sxy-timo', 'syy-timo');
end

figure;
plot(xpos(proby),ypos(proby), timosol_x(proby),timosol_y(proby));
legend('FDPM', 'Timo');

figure;
subplot(4,1, 1);
trisurf(tri, xpos,ypos, S11', 'FaceColor','interp'); view(0, 90);
subplot(4,1, 2);
trisurf(tri, xpos,ypos, S12', 'FaceColor','interp'); view(0, 90);
subplot(4,1, 3);
trisurf(tri, xpos,ypos, S21', 'FaceColor','interp'); view(0, 90);
subplot(4,1, 4);
trisurf(tri, xpos,ypos, S22', 'FaceColor','interp'); view(0, 90);

figure;
subplot(4,1, 1);
trisurf(tri, xpos,ypos, timosol_sxx, 'FaceColor','interp'); view(0, 90);
subplot(4,1, 2);
trisurf(tri, xpos,ypos, timosol_syy, 'FaceColor','interp'); view(0, 90);
subplot(4,1, 3);
trisurf(tri, xpos,ypos, timosol_sxy, 'FaceColor','interp'); view(0, 90);
% subplot(4,1, 4);
% trisurf(tri, xpos,ypos, S22', 'FaceColor','interp'); view(0, 90);

probx = find(nodeX == L0/2);
figure;
plot(nodeY(probx),S12(probx), 'x-', nodeY(probx),timosol_sxy(probx),'.-');
legend('sim-sxy','ana-sxy');



if 0
    Ktan = Ktan0;
    fdpmDriverShowLinElasticStiffEigen;
end





