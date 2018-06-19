
clear;


%% physical properties
par = struct;

% Young, Poisson
par.E = 100;
par.nu = 0.3;

% DP
par.c0 = 0.1;
[par.eta,par.etabar,par.xi] = materialDPSetAngle(20,10,'plain-strain');

% Cap
par.ptens = par.xi / par.eta * par.c0;
par.M = par.eta * sqrt(3);
par.beta = 0.25;
par.Hcap = 40;

p0 = 0;
p0 = -99999;
par.a0 = (par.ptens-p0) / (1 + par.beta);


%% geometry

R1 = 2;
H1 = 4;


% fun_u = @(x,y) 0.1 + 0.1*x + 0.2*y;
% fun_v = @(x,y) 0.05 + 0.15*x + 0.1*y;

fun_u = @(x,y) zeros(size(x));
fun_v = @(x,y) -0.0125*y;


%% fdpm setup

%
prob_type = 1; % plane-strain
% prob_type = 2; % axisymmetric
disp(['prob_type=',int2str(prob_type)]);

% particle number
nr1 = 10;
h0 = R1 / nr1;

% influence radius, if use p(2), must > 2
dilation = 2.1;
re = h0 * dilation;

% finite-increment-gradient stabilization
fig_stab_alpha = 0.0;
% fig_stab_alpha = 0.05;
% fig_stab_alpha = 0.5;
% fig_stab_alpha = 0.8;
% fig_stab_alpha = 1.0;

%
Estab = 0.0 * par.E;
% Estab = 1.0 * par.E;
Estab = 20.0 * par.E;
% Estab = 100.0 * par.E;
% Estab = 1000.0 * par.E;
% Estab = 10000.0 * par.E;
% Estab = 100000.0 * par.E;
Estab = Estab / h0^2;


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
    nh1 = round(H1/h0);
    for j = 1:nh1
    for i = 1:nr1
        xx = 0 + (i-0.5)*h0;
        yy = 0 + (j-0.5)*h0;
        
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
        % margin = 3;
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
        else
            % fixed_nodes_x(end+1,1) = numNodes;
            % fixed_disp_x(end+1,1) = 0;
            % fixed_nodes_y(end+1,1) = numNodes;
            % fixed_disp_y(end+1,1) = Htop * yy / (H1) + Hbot*(H1-yy)/H1;
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
    
    pause;
    % return;
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
[~,~,~,~,~,~,D,~] = materialDPC(par,prob_type, zeros(4,1), 0.0,0.0);

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

disp('Begin create initial stab matrix'); tic;
fdpmDriverCalcStab;
Kstab0 = Kstab;
Ktan0 = Ktan0 + Kstab0;
disp('End create initial stab matrix'); toc;

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
    
    return
end
































%% estimate tolerance
tol_rel = 1.0e-8;
tol_abs = 1.0e-8;
% tolerance based on loading force and displacement
if isempty(tracBCVals)
    tol1 = 0.0;
else
    tol1 = max(abs(tracBCVals));
end
% tolerance based on displacement
if isempty(dispBCDofs)
    tol2 = 0.0;
else
    tol2 = max(abs(Ktan0(dispBCDofs,dispBCDofs)*dispBCVals));
end
% tolearance for each iteration
tol = max(tol1,tol2) * tol_rel;
tol = max(tol, tol_abs);
disp(['tol = ', num2str(tol)]);

%% incremental loop

% total load step
maxStep = 1;
% maxStep = 10;
% maxStep = 20;

% max newton iteration
maxIter = 1;

% incremental loop
for step = 1:maxStep
    disp(['Load step = ', int2str(step)]);
    
    % external force
    fext = zeros(numDofs,1);
    fext(tracBCDofs) = (step/maxStep) * tracBCVals;
    
    % residual force
    fres = frct + fext - fint - fstab;
    
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
        
        % current position
        nodeDisp = reshape(uv, 2,numNodes);
        nodeCurr = nodeCoord + nodeDisp;
        
        % construct new stiffness
        Kassem = SpMatAssem(numDofs);
        fint = zeros(numDofs,1);
        
        for i = 1:numNodes
            nneigh = conn(i).numNeigh;
            ineigh = conn(i).neigh2;
            ineighdof = fdpmNodeDof(ineigh,'xy');
            
            % updated coords
            % ineighpos = nodeCoord(:,ineigh) + reshape(uv(ineighdof), 2,nneigh);
            ineighpos = nodeCurr(:,ineigh);
            
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
            
            %
            alpha_old = nodeAlpha_old(i);
            epbar_old = nodeEpbar_old(i);
            [dgama,dgamb,tau,epsE,epbar,alpha,D,mtype] = materialDPC(par,prob_type, epsEtr, epbar_old,alpha_old);
            
            % cauchy
            sigma = tau ./ detF;
            
            
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
            
            % save values
            nodeF(:,:,i) = Fi;
            nodeEpsE(:,i) = epsE;
            nodeEpbar(i) = epbar;
            nodeAlpha(i) = alpha;
            nodeSigma(:,i) = sigma;
            
        end % stiffness K matrix
        
        % create sparse K
        Ktan = Kassem.SpMatCreateAndClear();
        
        fdpmDriverCalcStab;
        Ktan = Ktan + Kstab;
        
        % update residual
        fres = frct + fext - fint - fstab;
        rnorm = norm(fres);
        disp(['|resid| = ', num2str(rnorm)]);
        if rnorm <= tol
            break;
        end
        
        
    end % Newton iteration
    
    % update states
    uv_old = uv;
    nodeEpsE_old = nodeEpsE;
    nodeF_old = nodeF;
    nodeEpbar_old = nodeEpbar;
    nodeAlpha_old = nodeAlpha;
    
    if 1
        nodeCurr = nodeCoord + reshape(uv, 2,numNodes);
        plot(nodePos(:,1),nodePos(:,2),'o', nodeCurr(1,:),nodeCurr(2,:),'x');
        legend('before','after');
        axis equal;
        axis([0 2 0 5]);
        title(['step=',int2str(step),'/',int2str(maxStep)]);
        
        drawnow;
    end
    
    
end % incremental loop


fdpmDriverWriteCsv;


return








