%%%testBarNecking


clear;



%% physical properties

% Young, Poisson
E0 = 206.9;
nu0 = 0.29;

% yield stress of vonMises
Fyield = 0.45;
% Fyield = 1.0;
% Fyield = 99999;

% bar size
R0 = 6.413e-3;
% Rmid = R0 * 1.0;
Rmid = R0 * 0.982;
% Rmid = R0 * 0.95;
% Rmid = R0 * 0.5;
L0 = 53.334e-3;
Lmid = L0 * 0.5;
Lref = 50.8e-3;
% Ldisp = 15.0e-3;
Ldisp = 5.0e-3;
% Ldisp = 2.0e-3;
% Ldisp = 1.0e-3;


%% method config

%
% prob_type = 1; % plane-strain
prob_type = 2; % axisymmetric
disp(['prob_type=',int2str(prob_type)]);

%
use_taylor_stab = 0;
% use_taylor_stab = 1;


%% fdpm setup

% particle number in radial direction 
nrad = 8;
h0 = R0/nrad;

nmid = round(Lmid / h0);
nlen = nmid * 2;

% influence radius, if use p(2), must > 2
% dilation = 1.5;
dilation = 1.8;
% dilation = 2.1;
% dilation = 2.5;
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
if 0
    testBarNeckingGenPart;
else
    testBarNeckingLoadMesh;
end

nodeX = nodeX + R0*0.00;

nodePos = [nodeX,nodeY];

% dof is (udisp,vdisp)
numDofs = numNodes * 2;

disp('End generate particles & BC');



if (1) % plot fdpm node mesh
    fdpmDriverPlotNodes;
end


%% build triangulation
if 1
    fdpmDriverBuildTri;
    % return;
end

%% build connnectivity
disp('Begin build connnectivity'); tic;
if 1
    % use_auto_re = 0;
    % fdpmDriverBuildConn;
    % fdpmDriverBuildConnIntegDomain;
    fdpmDriverBuildConnIntegBndry2;
    fdpmDriverBuildMoment;

    if prob_type==2 
        % nodeVol = nodeX .* nodeArea;
        nodeVol = 2*pi .* nodeVol;
    end
end
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
% 4x4, [11,22,12,33]
% D = E0/(1+nu0)/(1-2*nu0) * [...
% 1-nu0, nu0,   0,       nu0; 
% nu0,   1-nu0, 0,       nu0; 
% 0,     0,     0.5-nu0, 0; 
% nu0,   nu0,   0,       1-nu0];
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
Ktan = Ktan0;
disp('End create initial stiffness matrix'); toc;


%% estimate tolerance
tol_rel = 1.0e-9;
tol_abs = 1.0e-6;
% tol_abs = 1.0e-5;
% tolerance based on loading force and displacement
tol1 = max(abs(tracBCVals)); if isempty(tol1), tol1=0; end;
% tolerance based on displacement
tol2 = max(abs(Ktan0(dispBCDofs,dispBCDofs)*dispBCVals));
% tolearance for each iteration
tol = max(tol1,tol2) * tol_rel;
tol = max(tol, tol_abs);
disp(['tol = ', num2str(tol)]);

%% incremental loop

% total load step
maxStep = 500;
% maxStep = 200;
% max newton iteration
maxIter = 20;

% incremental loop
for step = 1:maxStep
    disp(['Load step = ', int2str(step)]);
    
    % external force
    fext = zeros(numDofs,1);
    % incremental loading
    fext(tracBCDofs) = (step/maxStep) * tracBCVals;
    
    % residual force
    fres = frct + fext - fint;
    
    % Newton-Raphson iteration
    for iter = 1:maxIter
        disp(['Newton iter = ', int2str(iter)]);
        
        if iter == 1
            % [dduv,dreact] = fdpmSolve(Ktan0,fres, dispBCDofs,(1/maxStep).*dispBCVals);
            [dduv,dreact] = fdpmSolve(Ktan,fres, dispBCDofs,(1/maxStep).*dispBCVals);
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
            % expand to full 3x3
            % Fi = fdpmFormFmat(Fi, prob_type, nodeCoord(1,i), ineighpos(1,end));
            if prob_type == 1 % plane strain, no deformation in Z-dir
                Fi(3,3) = 1;
            elseif prob_type == 2 % axisymmetric, hoop deform = x/X
                ri = nodeCoord(1,i);
                if ri == 0 % just on-axis
                    % Fi(3,3) = 1;
                    Fi(3,3) = dot(DN(1,:), ineighpos(1,:));
                else % off-axis, 
                    Fi(3,3) = ineighpos(1,end) / ri;
                end
            end
            
            % incremental defrom grad from previous step, F(n+1) = dF * F(n)
            dF = Fi / nodeF_old(:,:,i);
            
            % save deform grad
            detF = det(Fi);
            nodeF(:,:,i) = Fi;
            
            % current volume
            ivol = detF * nodeVol(i);
            
            % current d/dx
            dN = Fi(1:2,1:2)' \ DN;
            dNx = dN(1,:);
            dNy = dN(2,:);
            
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
            
            % material
            [D,tau,epsE,nodeFlag(i)] = materialVonMises(epsEtr, E0,nu0,Fyield, prob_type);
            
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
            
            if use_taylor_stab %
                ineighx = nodeX(ineigh);
                ineighy = nodeY(ineigh);
                % B = fdpmFormBmat(nneigh, conn(i).dNX,conn(i).dNY, prob_type, ineighx,ineighy);
                B = Bmat;
                Bx = fdpmFormBmat(nneigh, conn(i).dNXX, conn(i).dNXY, prob_type, ineighx,ineighy);
                By = fdpmFormBmat(nneigh, conn(i).dNXY, conn(i).dNYY, prob_type, ineighx,ineighy);
                Dmat = D;
                
                Kstab = Bx'*Dmat*Bx.*nodeMomXX(i) + By'*Dmat*By.*nodeMomYY(i);
                Kstab = Kstab + (Bx'*Dmat*By + By'*Dmat*Bx).*nodeMomXY(i);
                Kstab = Kstab + (B'*Dmat*Bx + Bx'*Dmat*B).*nodeMomX(i);
                Kstab = Kstab + (B'*Dmat*By + By'*Dmat*B).*nodeMomY(i);
                
                Ki = Ki + Kstab;
                
                Tstab = Bx.' * sigma(1:3) * nodeMomX(i) + By.' * sigma(1:3) * nodeMomY(i);
                Ti = Ti + Tstab;
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
    
    if 1
        dict = containers.Map;
        dict('sxx') = nodeSigma(1,:).';
        dict('syy') = nodeSigma(2,:).';
        dict('sxy') = nodeSigma(3,:).';
        dict('szz') = nodeSigma(4,:).';
        dict('disp') = reshape(uv,2,[]).';
        dict('epflag') = nodeFlag;
        
        
        vtk.writeTriMesh(sprintf('testBarNecking/output/hoge%04d.vtk',step), tri,nodeX,nodeY, dict);
    end

    
end % incremental loop



if 1
    dict = containers.Map;
    dict('sxx') = nodeSigma(1,:).';
    dict('syy') = nodeSigma(2,:).';
    dict('sxy') = nodeSigma(3,:).';
    dict('szz') = nodeSigma(4,:).';
    % dict('sxx_ana') = ana_sxx;
    % dict('syy_ana') = ana_syy;
    % dict('sxy_ana') = ana_sxy;
    dict('disp') = reshape(uv,2,[]).';
    % dict('disp_ana') = [ana_ux,ana_uy];
    dict('epflag') = nodeFlag;
    
    vtk.writeTriMesh('testBarNecking/hoge.vtk', tri,nodeX,nodeY, dict);
end




return








