
clear;


%% physical properties

% Young, Poisson
E0 = 100;
nu0 = 0.3;


%% geometry

R1 = 2;
H1 = 3;


fun_u = @(x,y) 0.1 + 0.1*x + 0.2*y;
fun_v = @(x,y) 0.05 + 0.15*x + 0.1*y;

% fun_u = @(x,y) zeros(size(x));
% fun_v = @(x,y) -0.0125*y;


%% fdpm setup

% coordinate system
prob_type = 1; % plane-strain
% prob_type = 2; % axisymmetric
disp(['prob_type=',int2str(prob_type)]);

% particle number
nr1 = 11;
h0 = R1 / (nr1-1);

% influence radius, if use p(2), must > 2
% dilation = 1.5;
dilation = 2.1;
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
    
    %
    nh1 = round(H1/h0) + 1;
    for j = 1:nh1
    for i = 1:nr1
        xx = 0 + (i-1)*h0;
        yy = 0 + (j-1)*h0;
        
        numNodes = numNodes + 1;
        nodeX(numNodes,1) = xx;
        nodeY(numNodes,1) = yy;
        
        vol = h0^2;
        if i==1 || i==nr1
            % vol = vol / 2;
        end
        if j==1 || j==nh1
            % vol = v / 2;
        end
        nodeVol(numNodes,1) = vol;
        
        margin = 1;
        if j <= margin % bottom
            fixed_nodes_y(end+1,1) = numNodes;
            fixed_nodes_x(end+1,1) = numNodes;
            fixed_disp_y(end+1,1) = fun_v(xx,yy);
            fixed_disp_x(end+1,1) = fun_u(xx,yy);
        elseif j >= nh1-margin+1 % top
            fixed_nodes_y(end+1,1) = numNodes;
            fixed_nodes_x(end+1,1) = numNodes;
            fixed_disp_y(end+1,1) = fun_v(xx,yy);
            fixed_disp_x(end+1,1) = fun_u(xx,yy);
        elseif i <= margin % left wall
            fixed_nodes_x(end+1,1) = numNodes;
            fixed_nodes_y(end+1,1) = numNodes;
            fixed_disp_y(end+1,1) = fun_v(xx,yy);
            fixed_disp_x(end+1,1) = fun_u(xx,yy);
        elseif i >= nr1-margin+1 % right wall
            fixed_nodes_x(end+1,1) = numNodes;
            fixed_nodes_y(end+1,1) = numNodes;
            fixed_disp_y(end+1,1) = fun_v(xx,yy);
            fixed_disp_x(end+1,1) = fun_u(xx,yy);
        end
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

if (1)
    % plot fdpm node mesh
    figure;
    plot(nodeX,nodeY,'o', ...
    nodeX(fixed_nodes_y),nodeY(fixed_nodes_y),'^', ...
    nodeX(fixed_nodes_x),nodeY(fixed_nodes_x),'>');
    legend('node','dispy','dispx');
    axis([0 2 0 5])
    axis('equal');
    
    % pause;
    % return;
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
fdpmDriverBuildConnIntegDomain;
disp('End build connnectivity'); toc;
% return;


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

nodeFlag = zeros(numNodes,1);

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

% D is the elastic tangent modulus
D = materialVonMises(zeros(4,1), E0,nu0, 1e9, prob_type); 


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

% create sparse K
Ktan0 = SpMatCreateAndClear(Kassem);

disp('End create initial stiffness matrix'); toc;

% disp('Begin create initial stab matrix'); tic;
% fdpmDriverCalcStab;
% Kstab0 = Kstab;
% Ktan0 = Ktan0 + Kstab0;
% disp('End create initial stab matrix'); toc;

if 1
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
    
    fdpmDriverWriteCsv;
    
end

if 1
    uv = reshape(uv, 2, []).';
    uana = fun_u(nodeX,nodeY);
    vana = fun_v(nodeX,nodeY);
    scale = 2.0;
    figure;
    triplot(tri, nodeX+scale.*uv(:,1),nodeY+scale.*uv(:,2));
    axis equal;
    hold on;
    plot(nodeX+scale.*uana, nodeY+scale.*vana,'rx');
    legend('sim','ana');
    hold off;
end


if 1
    Ktan = Ktan0;
    figure;
    neig = 10;
    [v,d] = eigs(Ktan, neig, 'sm');
    for ieig = 1:neig
        vv = v(:,ieig);
        vv = reshape(vv, 2, []).';
        triplot(tri, nodeX+scale.*vv(:,1), nodeY+scale.*vv(:,2));
        title(['eig(',int2str(ieig),')=',num2str(d(ieig,ieig))]);
        axis equal;
        pause;
    end
end








