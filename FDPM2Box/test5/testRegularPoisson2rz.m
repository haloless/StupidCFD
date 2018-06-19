%testRegularPoisson2rz
% Poisson 
%
% $$\nabla^2 u + f = 0$$ 
%

clear;



%% geometry
prob_type = 2; % RZ coord.


L = 1.0;

%% source function

srctype = 0;
if srctype == 0
    mm = 3;
    fun_fext = @(x,y) (4*pi^2 * x.^mm - mm^2 * x.^(mm-2)) .* sin(2*pi*y);
    fun_phi = @(x,y) x.^mm .* sin(2*pi*y);
% elseif srctype == 1
    % fun_fext = @(x,y) zeros(size(x));
    % fun_phi = @(x,y) 0.1 + 0.2*x + 0.3*y;
    % fun_phix = @(x,y) 0.2 * ones(size(x));
    % fun_phiy = @(x,y) 0.3 * ones(size(y));
% elseif srctype == 2
    % fun_fext = @(x,y) -8 .* ones(size(x));
    % fun_phi = @(x,y) x.^2 + 2*x.*y + 3 * y.^2;
    % fun_phix = @(x,y) 2*x + 2*y;
    % fun_phiy = @(x,y) 2*x + 6*y;
else
    error('Unknown problem');
end


%
% meshtype = 0; % regular
meshtype = 1; % perturb
if meshtype == 1
    rng(0); % reset random seed
end

%% fdpm setup

% nlen = 8;
% nlen = 11;
nlen = 16;
% nlen = 21;
% nlen = 26;
% nlen = 31;
% nlen = 36;
% nlen = 46;
% nlen = 51;
nlen1 = nlen - 1;
h0 = L / nlen1;


% influence radius, if use p(2), must > 2
% dilation = 1.1;
% dilation = 1.5;
% dilation = 1.6;
dilation = 1.8;
% dilation = 2.1;
% dilation = 3.1;
re = h0 * dilation;


%% generate points and boundary conditions

% generation routine should setup the following data
numNodes = 0;
% particle states
nodeX = [];
nodeY = [];
nodeVol = [];
% BC Dirichlet
dirBCDofs = [];
dirBCVals = [];

disp('Begin generate particles & BC');
if 1
    for j = 1:nlen
    for i = 1:nlen
        xx = (i-1) * h0;
        yy = (j-1) * h0;
        
        if meshtype == 1 % perturb internal points, not boundary points
            if i~=1 && i~=nlen && j~=1 && j~=nlen
                rr = 0.4;
                xx = xx + (rand(1)-0.5)*h0*rr;
                yy = yy + (rand(1)-0.5)*h0*rr;
            end
        end
        
        nn = numNodes + 1;
        numNodes = nn;
        nodeX(nn,1) = xx;
        nodeY(nn,1) = yy;
        
        vv = h0^2;
        if i==1 || i==nlen
            vv = vv / 2;
        end
        if j==1 || j==nlen
            vv = vv / 2;
        end
        nodeVol(nn,1) = vv;
        
        % BC
        if i==1 || i==nlen || j==1 || j==nlen
            dirBCDofs(end+1,1) = nn;
            dirBCVals(end+1,1) = fun_phi(xx,yy);
        end
    end
    end
end
disp('End generate particles & BC');

nodePos = [nodeX,nodeY];

% 
numDofs = numNodes * 1;

disp('Begin create boundary');

disp('End create boundary');
nbndry = 0;

if (1)
    % plot fdpm node mesh
    figure;
    plot(nodeX,nodeY,'o', ...
    nodeX(dirBCDofs),nodeY(dirBCDofs),'x');
    legend('node','dir');
    % axis([0 L 0 L]);
    axis('equal');
    
    % pause;
    % return;
end

%
if (1)
    fdpmDriverBuildTri;
    nodeVol = nodeArea;
    % pause;
end

%% build connnectivity

disp('Begin build connnectivity'); tic;
if 0
    fdpmDriverBuildConnRegular;
else
    % fdpmDriverBuildConnIntegDomain;
    % fdpmDriverBuildConnIntegBndry;
    fdpmDriverBuildConnIntegBndry2;
    fdpmDriverBuildMoment;
end
disp('End build connnectivity'); toc;
% return

%% allocate buffers

fext = zeros(numDofs,1);
fint = zeros(numDofs,1);

%% create initial stiffness matrix
disp('Begin create initial stiffness matrix'); tic;

% assembler
Kassem = SpMatAssem(numDofs);

for i = 1:numNodes
    nneigh = conn(i).numNeigh;
    ineigh = conn(i).neigh2;
    ivol = nodeVol(i);
    
    dN = [conn(i).dNX; conn(i).dNY];
    Ke = dN' * dN * ivol;
    
    SpMatAssemBlockWithDof(Kassem, Ke,ineigh,ineigh);
end

if 1
    fdpmDriverPoissonStabHessian;
end

% create sparse K
Ktan0 = SpMatCreateAndClear(Kassem);

disp('End create initial stiffness matrix'); toc;



disp('Begin solve Poisson');

% RHS vector
fext = fun_fext(nodeX,nodeY);
fext = fext .* nodeVol;

% solve
[phi] = fdpmSolve(Ktan0, fext, dirBCDofs,dirBCVals);

disp('End solve Poisson');


if 1
    figure;
    % xx = reshape(nodeX, nlen,nlen);
    % yy = reshape(nodeY, nlen,nlen);
    % ss = reshape(phi, nlen,nlen);
    % surf(xx,yy,ss, 'FaceColor','interp');
    trisurf(tri, nodeX,nodeY,phi, 'FaceColor','interp');
    % zlim([-0.06 0.06]);
    title('solution');
    
    figure;
    % zz = fun_phi(xx,yy);
    zz = fun_phi(nodeX,nodeY);
    % surf(xx,yy, zz, 'FaceColor','interp');
    trisurf(tri, nodeX,nodeY, zz, 'FaceColor','interp');
    % zlim([-0.06 0.06])
    title('analytical');
end


if 0
    Ktan = Ktan0;
    fdpmDriverShowStiffEigen;
end


if 1
    fdpmDriverShowIntegError;
    vsum0 = v_sum;
end


if 0
    Ktan = Ktan2;
    fdpmDriverShowStiffEigen;
end



return






