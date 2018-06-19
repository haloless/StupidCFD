% Solve thick wall tube problem
% Use Lagrange Multiplier for BC


clear;



%% physical properties

% Young, Poisson
E0 = 210;
nu0 = 0.3;

% internal pressure
P0 = 1.0;
% P0 = 50.0;

% inner/outer radius
Ra = 0.1;
Rb = 0.2;
% length, not used for plane-strain
H = 0.1;


%% fdpm setup

%
prob_type = 1; % plane strain
% prob_type = 2; % axisymmetric
% assert(prob_type == 1);


% particle number in radial direction 
nrad = 20;
% nrad = 40;
% nrad = 80;
% nrad = 100;
h0 = (Rb-Ra)/nrad;

% influence radius, if use p(2), must > 2
dilation = 2.1;
% dilation = 3.1;
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
if prob_type == 1
    % plane strain
    % generate a circular slice of the tube
    
    % used to find most inner particles for apply pressure
    loaded_nodes = [];
    loaded_force = [];
    
    for irad = 1:nrad
        if 0
            % 1/4 cylinder
            rr = Ra + (irad-0.5)*h0;
            ll = pi/2 * rr;
            nn = round(ll/h0);
            vv = pi/4 * ((rr+0.5*h0)^2-(rr-0.5*h0)^2) / (nn);
            
            for ii = 1:nn
                tt = pi/2 / (nn) * (ii-0.5);
                
                numNodes = numNodes + 1;
                nodeX(numNodes,1) = rr * cos(tt);
                nodeY(numNodes,1) = rr * sin(tt);
                nodeVol(numNodes,1) = vv;
                
                if irad == 1
                    pp = P0*(pi/2*Ra)/(nn) * [cos(tt),sin(tt)];
                    
                    loaded_nodes(end+1,1) = numNodes;
                    loaded_force(1:2,end+1) = pp;
                end
            end
        else
            % full cylinder
            rr = Ra + (irad-0.5)*h0;
            ll = pi*2 * rr;
            nn = round(ll/h0);
            vv = pi * ((rr+0.5*h0)^2-(rr-0.5*h0)^2) / (nn);
            
            for ii = 1:nn
                tt = pi*2 / (nn) * (ii-0.5);
                
                numNodes = numNodes + 1;
                nodeX(numNodes,1) = rr * cos(tt);
                nodeY(numNodes,1) = rr * sin(tt);
                nodeVol(numNodes,1) = vv;
                
                if irad == 1
                    pp = P0*(pi*2*Ra)/(nn) * [cos(tt),sin(tt)];
                    
                    loaded_nodes(end+1,1) = numNodes;
                    loaded_force(1:2,end+1) = pp;
                end
            end
        end
    end
    
    % this setup does not use fixed dof
    
    % setup Loading
    tracBCDofs = fdpmNodeDof(loaded_nodes,'xy');
    tracBCVals = reshape(loaded_force, [],1);
    
elseif prob_type == 2
    % axisymmetric
    % simply generate a rectangular block
    
    % particle number in height direction
    nlen = round(H/h0);
    
    % used to find most inner particles for apply pressure
    loaded_nodes = [];
    loaded_force = [];
    
    for irad = 1:nrad
    for jlen = 1:nlen
        if 1
            rr = (irad-0.5)*h0 + Ra;
            zz = (jlen-0.5)*h0 + 0;
            
            numNodes = numNodes + 1;
            
            nodeX(numNodes,1) = rr;
            nodeY(numNodes,1) = zz;
            
            nodeVol(numNodes,1) = 2*pi*rr * h0^2;
            
            if irad == 1
                pp = P0*(2*pi*Ra)*h0;
                if jlen==1 || jlen==nlen
                    % pp = pp * 0.5;
                end
                
                loaded_nodes(end+1,1) = numNodes;
                loaded_force(1:2,end+1) = [pp,0];
            end
        end
    end
    end
    
    % do not need Dirichlet BC
    
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

%% create constraints
disp('Begin generate constraints'); tic;
numCons = 0;
consX = [];
consY = [];
consU = [];
consV = [];
if prob_type == 1
    % constrain on x-axis
    for irad = 0:nrad
        rr = Ra + (irad)*h0;
        % tt = 0;
        tt = pi/4;
        numCons = numCons + 1;
        consX(numCons,1) = rr * cos(tt);
        consY(numCons,1) = rr * sin(tt);
        consU(numCons,1) = sin(tt);
        consV(numCons,1) = -cos(tt);
    end
    % constrain on y-axis
    for irad = 0:nrad
        rr = Ra + (irad)*h0;
        % tt = pi/2;
        tt = pi/4*3;
        numCons = numCons + 1;
        consX(numCons,1) = rr * cos(tt);
        consY(numCons,1) = rr * sin(tt);
        consU(numCons,1) = sin(tt);
        consV(numCons,1) = -cos(tt);
    end
elseif prob_type == 2
    % two constrained on x-axis for bottom and top
    for irad = 0:nrad
        numCons = numCons + 1;
        consX(numCons,1) = Ra + irad*h0;
        consY(numCons,1) = 0;
        consU(numCons,1) = 0;
        consV(numCons,1) = 1;
    end
    for irad = 0:nrad
        numCons = numCons + 1;
        consX(numCons,1) = Ra + irad*h0;
        consY(numCons,1) = H;
        consU(numCons,1) = 0;
        consV(numCons,1) = 1;
    end
end
consPos = [consX,consY];
disp('End generate constraints'); toc;



if (1)
    % plot fdpm node mesh
    figure;
    % nodeX(fixed_nodes_y),nodeY(fixed_nodes_y),'^', 
    % nodeX(fixed_nodes_x),nodeY(fixed_nodes_x),'>', 
    plot(nodeX,nodeY,'o', ...
    nodeX(loaded_nodes),nodeY(loaded_nodes),'s', ...
    consX,consY,'.');
    legend('node','trac','cons');
    axis('equal');
    hold on;
    
    % plot load
    quiver(nodeX(loaded_nodes),nodeY(loaded_nodes),loaded_force(1,:)',loaded_force(2,:)');
    
    % plot constraint
    quiver(consX,consY,consU,consV);
    
    hold off;
    
    % return;
end


%% build connnectivity

disp('Begin build connnectivity'); tic;
fdpmDriverBuildConn;
disp('End build connnectivity'); toc;

%% build constraint
disp('Begin build constraint'); tic;
clear cons;
for k = 1:numCons
    [neigh,rs,rx,ry,cutoff] = fdpmNeighborhood(consPos(k,:),nodePos,re);
    
    % weight function
    w = fdpmWeightFunc(rs(neigh),cutoff);
    
    % shape function
    Nshape = fdpmShapeMLS(rx(neigh),ry(neigh),w);
    
    % save 
    cons(k).numNeigh = numel(neigh);
    cons(k).neigh = neigh;
    cons(k).Nshape = Nshape;
end
disp('End build constraint'); toc;

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
% 3x3, [11,22,12]
% D = E0/(1+nu0)/(1-2*nu0) * [1-nu0,nu0,0;nu0,1-nu0,0;0,0,0.5-nu0];
% 4x4, [11,22,12,33]
D = E0/(1+nu0)/(1-2*nu0) * [1-nu0,nu0,0,nu0; nu0,1-nu0,0,nu0; 0,0,0.5-nu0,0; nu0,nu0,0,1-nu0];

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

% create Lagrange multipler constraint matrix
disp('Begin create constraint matrix'); tic;
Kcons = sparse(numDofs,numCons);
for k = 1:numCons
    kneigh = cons(k).neigh;
    
    kneighdof = fdpmNodeDof(kneigh,'x');
    Kcons(kneighdof,k) = Kcons(kneighdof,k) - consU(k).*cons(k).Nshape';
    
    kneighdof = fdpmNodeDof(kneigh,'y');
    Kcons(kneighdof,k) = Kcons(kneighdof,k) - consV(k).*cons(k).Nshape';
end
disp('End create constraint matrix'); toc;


if 1
    % solve infinitesimal problem and goodbye
    fext = zeros(numDofs,1);
    fext(tracBCDofs) = tracBCVals;
    
    frhs = fext - fint;
    
    % sol = fdpmSolve(Ktan0,frhs, dispBCDofs,dispBCVals);
    
    Kall = [Ktan0, Kcons; Kcons', sparse(numCons,numCons)];
    ball = [frhs; zeros(numCons,1)];
    sol = Kall \ ball;
    sol = sol(1:numDofs);
    
    uv = reshape(sol, 2,numNodes);
    
    % update particle coordinates
    nodeCoord = nodeCoord + uv * 10;
    % xpos = reshape(nodeCoord(1,:), nx,ny);
    % ypos = reshape(nodeCoord(2,:), nx,ny);
    
    % plot particles after displacement
    figure;
    plot(nodePos(:,1),nodePos(:,2),'o', nodeCoord(1,:),nodeCoord(2,:),'x');
    legend('before','after');
    axis equal;
    
    
    % plot radial displacement
    figure;
    if prob_type == 1
        % plane-strain config
        % average radial disp for each layer
        rr2 = sqrt(nodeX.^2+nodeY.^2)';
        ur2 = (uv(1,:).*nodeX' + uv(2,:).*nodeY')./rr2;
        for irad = 1:nrad
            rr(irad) = Ra + (irad-0.5)*h0;
            ur(irad) = mean(ur2(abs(rr2-rr(irad))<1.0e-6));
        end
    elseif prob_type == 2
        % axisymmetric config
        rr2 = nodeX.';
        ur2 = uv(1,:);
        for irad = 1:nrad
            rr(irad) = Ra + (irad-0.5)*h0;
            ur(irad) = mean(ur2(abs(rr2-rr(irad))<1.0e-6));
        end
    end
    
    % analytical plane-strain solution
    c1 = (1+nu0)/E0*(1-2*nu0)*P0/(Rb^2/Ra^2-1);
    c2 = (1+nu0)/E0*Rb^2*P0/(Rb^2/Ra^2-1);
    rs = linspace(Ra,Rb,100);
    us = c1 .* rs + c2 ./ rs;
    plot(rr,ur,'x',rs,us,'-');
    
    return;
end

%% estimate tolerance
tol_rel = 1.0e-9;
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
maxStep = 1;
% maxStep = 2;
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
            
            % trial e
            etr = 0.5 * logm(BeTr);
            epsEtr = [ etr(1); etr(5); etr(2)*2; etr(9) ];
            
            % TODO material routine
            epsE = epsEtr;
            if prob_type == 1
                assert(abs(epsE(4))<1.0e-9); % out-of-plane strain should be essentially zero
            end
            nodeEpsE(:,i) = epsE;
            
            % kirchhoff
            tau = D * epsE;
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
    
    nodeCoord = nodeCoord + uv;
    
    figure;
    plot(nodePos(:,1),nodePos(:,2),'o', nodeCoord(1,:),nodeCoord(2,:),'x');
    legend('before','after');
    axis equal;
    
    
    figure;
    if prob_type == 1
        rr2 = sqrt(nodeX.^2+nodeY.^2)';
        ur2 = (uv(1,:).*nodeX' + uv(2,:).*nodeY')./rr2;
    elseif prob_type == 2
        rr2 = nodeX.';
        ur2 = uv(1,:);
    end
    for irad = 1:nrad
        rr(irad) = Ra + (irad-0.5)*h0;
        ur(irad) = mean(ur2(abs(rr2-rr(irad))<1.0e-6));
    end
    
    % analytical plane-strain solution
    c1 = (1+nu0)/E0*(1-2*nu0)*P0/(Rb^2/Ra^2-1);
    c2 = (1+nu0)/E0*Rb^2*P0/(Rb^2/Ra^2-1);
    rs = linspace(Ra,Rb,100);
    us = c1 .* rs + c2 ./ rs;
    plot(rr,ur,'x',rs,us,'-');
end


return








