%testPlateHole: elastic infinite plate with a hole, Timoshenko




clear;

%% physical property

prob_type = 1; % plane-strain

E = 100;
nu = 0.3;

R = 1.0; % hole radius
Tx = 1.0;

% analytical solution
ana = refsol_PlateWithHole(E, nu, R, Tx);


%% generate points and boundary conditions

% generation routine should setup the following data
numNodes = 0;
% particle states
nodeX = [];
nodeY = [];
% BC
dispBCDofs = [];
dispBCVals = [];
tracBCDofs = [];
tracBCVals = [];

disp('Begin generate particles & BC');

% load gmsh
[groups,nodes] = gmsh.loadMeshFile('plate.msh');

% node coordinates
numNodes = size(nodes,1);
nodeX = nodes(:,1);
nodeY = nodes(:,2);
nodePos = [nodeX,nodeY];

numDofs = numNodes * 2;



% boundary condition
loaded_nodes = [];
loaded_force = [];
fixed_nodes_x = [];
fixed_disp_x = [];
fixed_nodes_y = [];
fixed_disp_y = [];

if 0 % set all displacement
    % 
    ids = groups('left').nodeIds; num = length(ids);
    fixed_nodes_x(end+1:end+num,1) = ids;
    fixed_disp_x(end+1:end+num,1) = zeros(num,1);
    % 
    ids = groups('bottom').nodeIds; num = length(ids);
    fixed_nodes_y(end+1:end+num,1) = ids;
    fixed_disp_y(end+1:end+num,1) = zeros(num,1);
    %
    ids = unique([groups('right').nodeIds; groups('top').nodeIds]);
    num = length(ids);
    fixed_nodes_x(end+1:end+num,1) = ids;
    fixed_disp_x(end+1:end+num,1) = ana.ux(nodeX(ids),nodeY(ids));
    fixed_nodes_y(end+1:end+num,1) = ids;
    fixed_disp_y(end+1:end+num,1) = ana.uy(nodeX(ids),nodeY(ids));
else % set displacement and traction
    % 
    ids = groups('left').nodeIds; num = length(ids);
    fixed_nodes_x(end+1:end+num,1) = ids;
    fixed_disp_x(end+1:end+num,1) = zeros(num,1);
    % 
    ids = groups('bottom').nodeIds; num = length(ids);
    fixed_nodes_y(end+1:end+num,1) = ids;
    fixed_disp_y(end+1:end+num,1) = zeros(num,1);
    
    loaded_nodes = zeros(numNodes,1);
    loaded_force = zeros(2,numNodes);
    
    %
    ids = groups('right').conn;
    for k = 1:size(ids,1)
        k1 = ids(k,1); k2 = ids(k,2);
        x1 = nodeX(k1); y1 = nodeY(k1);
        x2 = nodeX(k2); y2 = nodeY(k2);
        len = norm([x2-x1; y2-y1]);
        
        loaded_nodes(k1) = 1;
        loaded_force(:,k1) = loaded_force(:,k1) + [ana.sxx(x1,y1); ana.sxy(x1,y1)] .* (len/2);
        
        loaded_nodes(k2) = 1;
        loaded_force(:,k2) = loaded_force(:,k2) + [ana.sxx(x2,y2); ana.sxy(x2,y2)] .* (len/2);
    end
    %
    ids = groups('top').conn;
    for k = 1:size(ids,1)
        k1 = ids(k,1); k2 = ids(k,2);
        x1 = nodeX(k1); y1 = nodeY(k1);
        x2 = nodeX(k2); y2 = nodeY(k2);
        len = norm([x2-x1; y2-y1]);
        
        loaded_nodes(k1) = 1;
        loaded_force(:,k1) = loaded_force(:,k1) + [ana.sxy(x1,y1); ana.syy(x1,y1)] .* (len/2);
        
        loaded_nodes(k2) = 1;
        loaded_force(:,k2) = loaded_force(:,k2) + [ana.sxy(x2,y2); ana.syy(x2,y2)] .* (len/2);
    end
    
    loaded_nodes = find(loaded_nodes == 1);
    loaded_force = loaded_force(:,loaded_nodes);
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


disp('End generate particles & BC');

if (1) % plot fdpm node mesh
    fdpmDriverPlotNodes;
end




%% fdpm setup

h0 = 0.2; % TODO estimate from mesh file

% influence radius, if use p(2), must > 2
% dilation = 1.1;
% dilation = 1.5;
% dilation = 1.6;
dilation = 1.8;
% dilation = 2.1;
% dilation = 3.1;
re = h0 * dilation;



nodeVol = [];
%
if (1)
    edgePair = [];
    edgePair = [edgePair; groups(1).conn];
    edgePair = [edgePair; groups(2).conn];
    edgePair = [edgePair; groups(3).conn];
    edgePair = [edgePair; groups(4).conn];
    edgePair = [edgePair; groups(5).conn];
    numEdges = size(edgePair,1);
    fdpmDriverBuildTri;
    nodeVol = nodeArea;
    % pause;
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

% elastic 
D = materialVonMises(zeros(4,1), E,nu, -1, prob_type); 
if prob_type == 1
    Dmat = D(1:3,1:3);
elseif prob_type == 2
    Dmat = D;
end


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

if 1
    nodeFlag2 = nodeFlag;
    % for i = 1:numTri
        % if any(nodeFlag(tri(i,:)))
            % nodeFlag2(tri(i,:)) = 1;
        % end
    % end
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

[ana_ux,ana_uy] = ana.disp_xy(nodeX, nodeY);
[ana_sxx,ana_syy,ana_sxy] = ana.stress(nodeX, nodeY);


if 1
    figure;
    trisurf(tri, nodeX,nodeY, nodeSigma(1,:)./Tx, 'FaceColor','interp');
    view(0, 90); axis equal; xlabel('x'); ylabel('y'); title('sim-sxx');
end

if 1
    [~,yprobe] = sort(nodeY(groups('left').nodeIds));
    yprobe = groups('left').nodeIds(yprobe);
    figure;
    plot(nodeY(yprobe),nodeSigma(1,yprobe),'x-', nodeY(yprobe),ana_sxx(yprobe),'.-');
    
end



if 1
    % figure;
    % trisurf(tri, nodeX,nodeY, ana_ux, 'FaceColor','interp');
    % view(0, 90); axis equal; xlabel('x'); ylabel('y'); title('ux');
    % figure;
    % trisurf(tri, nodeX,nodeY, ana_uy, 'FaceColor','interp');
    % view(0, 90); axis equal; xlabel('x'); ylabel('y'); title('uy');
    figure;
    trisurf(tri, nodeX,nodeY, ana_sxx./Tx, 'FaceColor','interp');
    view(0, 90); axis equal; xlabel('x'); ylabel('y'); title('ana-sxx');
    
end

if 1
    dict = containers.Map;
    dict('sxx') = nodeSigma(1,:).';
    dict('syy') = nodeSigma(2,:).';
    dict('sxy') = nodeSigma(3,:).';
    dict('sxx_ana') = ana_sxx;
    dict('syy_ana') = ana_syy;
    dict('sxy_ana') = ana_sxy;
    dict('disp') = nodeDisp.';
    dict('disp_ana') = [ana_ux,ana_uy];
    
    vtk.writeTriMesh('testPlateHole/hoge00.vtk', tri,nodeX,nodeY, dict);
end




%% 





