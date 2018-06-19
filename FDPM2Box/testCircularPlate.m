% Plastic FEM book, Sec. 7.5.3, P280

clear;


%% physical properties

% Young, Poisson
E0 = 1.0e7;
nu0 = 0.24;
% yield stress of vonMises
Fyield = 16000;

% load pressure
% P0 = 100;
% P0 = 200;
% P0 = 250;
% P0 = 260;
% P0 = 270;
% P0 = 275;
% P0 = 280;
P0 = 285; % almost critical
% P0 = 300;


% inner/outer radius
R = 10;
% length
H = 1.0;


%% fdpm setup

%
prob_type = 2; % this problem is axisymmetric
prob_type

% particle number in height direction 
nlen = 16;
h0 = H / nlen;
nrad = round(R/h0);

% influence radius, if use p(2), must > 2
dilation = 2.1;
re = h0 * dilation;

% finite-increment-gradient stabilization
% fig_stab_alpha = 0.0;
fig_stab_alpha = 0.5;


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
if prob_type == 2
    % axisymmetric
    % simply generate a rectangular block
    
    % used to find most inner particles for apply pressure
    loaded_nodes = [];
    loaded_force = [];
    % used to find lower/upper particles to constraint on x-axis
    fixed_nodes_x = [];
    fixed_nodes_y = [];
    
    for jlen = 1:nlen
    for irad = 1:nrad
        rr = 0 + (irad-0.5)*h0;
        zz = 0 + (jlen-0.5)*h0;
        
        numNodes = numNodes + 1;
        
        nodeX(numNodes,1) = rr;
        nodeY(numNodes,1) = zz;
        
        nodeVol(numNodes,1) = 2*pi*rr * h0^2;
        
        if jlen == nlen
            pp = P0 * 2*pi*rr*h0;
            loaded_nodes(end+1,1) = numNodes;
            loaded_force(1:2,end+1) = [0,-pp];
        end
        
        if jlen==1 && irad==nrad
            fixed_nodes_y(end+1,1) = numNodes;
            fixed_nodes_x(end+1,1) = numNodes;
        elseif irad==1
            % fixed_nodes_x(end+1,1) = numNodes;
        end
    end
    end
    
    % setup Dirichlet BC
    tmp = fdpmNodeDof(fixed_nodes_y,'y');
    dispBCDofs = [dispBCDofs; tmp];
    dispBCVals = [dispBCVals; zeros(size(tmp))];
    tmp = fdpmNodeDof(fixed_nodes_x,'x');
    dispBCDofs = [dispBCDofs; tmp];
    dispBCVals = [dispBCVals; zeros(size(tmp))];
    
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


if (0)
    % plot fdpm node mesh
    figure;
    % fdpmPlotNodes;
    plot(nodeX,nodeY,'o', ...
    nodeX(fixed_nodes_y),nodeY(fixed_nodes_y),'^', ...
    nodeX(fixed_nodes_x),nodeY(fixed_nodes_x),'>', ...
    nodeX(loaded_nodes),nodeY(loaded_nodes),'s');
    legend('node','dispy','dispx','trac');
    axis('equal');
    hold on;
    quiver(nodeX(loaded_nodes),nodeY(loaded_nodes),loaded_force(1,:)',loaded_force(2,:)');
    hold off;
    
    return;
end


%% build connnectivity

disp('Begin build connnectivity'); tic;
fdpmDriverBuildConn;
disp('End build connnectivity'); toc;


%% allocate buffers

% [fx,fy]
fext = zeros(numDofs,1);
fint = zeros(numDofs,1);
frct = zeros(numDofs,1);

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

nodeFlag = zeros(numNodes,1);

% initialize
for i = 1:numNodes
    nodeF(:,:,i) = eye(3);
end

% save old values
uv_old = uv;
nodeEpsE_old = nodeEpsE;
nodeF_old = nodeF;

%% create initial stiffness matrix
disp('Begin create initial stiffness matrix'); tic;

% assembler
Kassem = SpMatAssem(numDofs);

% D is the elastic tangent modulus
D = materialVonMises(zeros(4,1), E0,nu0,Fyield, prob_type);

for i = 1:numNodes
    nneigh = conn(i).numNeigh;
    ineigh = conn(i).neigh2;
    ineighdof = fdpmNodeDof(ineigh,'xy');
    ivol = nodeVol(i);
    
    B = fdpmFormBmat(nneigh, conn(i).dNX,conn(i).dNY, prob_type,nodeX(ineigh),nodeY(ineigh));
    
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


%% estimate tolerance
tol_rel = 1.0e-6;
tol_abs = 1.0e-6;
% tolerance based on loading force and displacement
tol1 = max(abs(tracBCVals));
% tolerance based on displacement
tol2 = max(abs(Ktan0(dispBCDofs,dispBCDofs)*dispBCVals));
% tolearance for each iteration
tol = max(tol1,tol2) * tol_rel;
tol = max(tol, tol_abs);
disp(['tol = ', num2str(tol)]);

%% incremental loop

% total load step
% maxStep = 5;
maxStep = 10;
% max newton iteration
maxIter = 40;

% incremental loop
for step = 1:maxStep
    disp(['Load step = ', int2str(step)]);
    
    % external force
    fext = zeros(numDofs,1);
    fext(tracBCDofs) = (step/maxStep) * tracBCVals;
    
    % residual force
    fres = frct + fext - fint;
    
    % Newton-Raphson iteration
    for iter = 1:maxIter
        disp(['Newton iter = ', int2str(iter)]);
        
        if iter == 1
            [dduv,dreact] = fdpmSolve(Ktan0,fres, dispBCDofs,(1/maxStep).*dispBCVals);
        else
            [dduv,dreact] = fdpmSolve(Ktan, fres, dispBCDofs,[]);
        end
        
        % update
        uv = uv + dduv;
        frct = frct + dreact;
        duv = uv - uv_old;
        
        % construct new stiffness
        Kassem = SpMatAssem(numDofs);
        fint = zeros(numDofs,1);
        
        for i = 1:numNodes
            nneigh = conn(i).numNeigh;
            ineigh = conn(i).neigh2;
            ineighdof = fdpmNodeDof(ineigh,'xy');
            
            % updated coords
            ineighpos = nodeCoord(:,ineigh) + reshape(uv(ineighdof), 2,nneigh);
            
            % d/dX
            DN = [conn(i).dNX; conn(i).dNY];
            
            % total deform grad in-plane
            Fi = ineighpos * DN';
            
            % d/dx
            dN = Fi' \ DN;
            dNx = dN(1,:);
            dNy = dN(2,:);
            
            % expand to full 3x3
            if prob_type == 1
                Fi(3,3) = 1.0;
            elseif prob_type == 2
                Fi(3,3) = ineighpos(1,end) / nodeCoord(1,i);
            end
            detF = det(Fi);
            nodeF(:,:,i) = Fi;
            
            % current volume
            ivol = detF * nodeVol(i);
            
            % incremental defrom grad from previous step, F(n+1) = dF * F(n)
            dF = Fi / nodeF_old(:,:,i);
            
            % previous elastic log-strain
            e = nodeEpsE_old(:,i);
            eps = [e(1),e(3)/2,0; e(3)/2,e(2),0; 0,0,e(4)];
            
            % trial Be
            Be = expm(2*eps);
            BeTr = dF * Be * dF.';
            
            % trial log-strain e
            % etr = 0.5 * logm(BeTr);
            etr = 0.5 * fdpmLogm(BeTr);
            epsEtr = [ etr(1); etr(5); etr(2)*2; etr(9) ];
            
            if 0
                % TODO material routine
                epsE = epsEtr;
                if prob_type == 1
                    assert(abs(epsE(4))<1.0e-9); % out-of-plane strain should be essentially zero
                end
                
                tau = D * epsE;
            else
                [D,tau,epsE,nodeFlag(i)] = materialVonMises(epsEtr, E0,nu0,Fyield, prob_type);
            end
            
            nodeEpsE(:,i) = epsE;
            
            % cauchy
            sigma = tau ./ detF;
            nodeSigma(:,i) = sigma;
            
            % spatial tangent modulus
            amat = fdpmFormAmat(BeTr, sigma, D, detF, prob_type);
            
            % gradient operator
            Bmat = fdpmFormBmat(nneigh, dNx,dNy, prob_type,ineighpos(1,:),ineighpos(2,:));
            Gmat = fdpmFormGmat(nneigh, dNx,dNy, prob_type,ineighpos(1,:),ineighpos(2,:));
            
            % point stiffness
            Ki = Gmat.' * amat * Gmat * ivol;
            
            % point force
            if prob_type == 1
                Ti = Bmat.' * sigma(1:3) * ivol;
            elseif prob_type == 2
                Ti = Bmat.' * sigma * ivol;
            end
            
            % assemble
            Kassem.SpMatAssemBlockWithDof(Ki,ineighdof,ineighdof);
            fint(ineighdof) = fint(ineighdof) + Ti;
        end
        
        
        % update residual
        fres = frct + fext - fint;
        rnorm = norm(fres);
        disp(['|resid| = ', num2str(rnorm)]);
        if rnorm <= tol
            break;
        end
        
        % create sparse K
        Ktan = Kassem.SpMatCreateAndClear();
        
    end % Newton iteration
    
    % update states
    uv_old = uv;
    nodeEpsE_old = nodeEpsE;
    nodeF_old = nodeF;
    
end % incremental loop


if 1
    uv = reshape(uv, 2,numNodes);
    nodeCurr = nodeCoord + uv;
    
    % current config
    figure;
    plot(nodePos(:,1),nodePos(:,2),'o', nodeCurr(1,:),nodeCurr(2,:),'x');
    legend('before','after');
    axis equal;
    
    
    figure;
    epmask = nodeFlag==1;
    plot(nodeCurr(1,epmask),nodeCurr(2,epmask),'x', nodeCurr(1,~epmask),nodeCurr(2,~epmask),'o')
    axis equal;
    legend('plastic','elastic');
    
    figure;
    for irad = 1:nrad
        rr = 0 + (irad-0.5)*h0;
        ii = find(abs(nodePos(:,1)-rr)<1.0e-6);
        vv = mean(uv(2,ii));
        % vv = min(uv(2,ii));
        
        rr2(irad) = rr;
        vv2(irad) = vv;
    end
    plot(rr2./R,vv2./H,'x-');
    xlabel('r/R'); ylabel('v/H');
    axis([0,1, -0.6,0.0]);
end

return








