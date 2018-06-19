%testRegularPoisson2
% Poisson 
%
% $$\nabla^2 u + f = 0$$ 
%

clear;



%% geometry
L = 1.0;

%% source function

srctype = 0;
% srctype = 1;
% srctype = 2;
if srctype == 0
    mm = 1; nn = 1;
    % mm = 2; nn = 3;
    % f=sin(m*pi*x) * sin(m*pi*y)
    fun_fext = @(x,y) (mm^2+nn^2)/2 * sin(mm*pi*x) .* sin(nn*pi*y);
    fun_phi = @(x,y) 1/(2*pi^2) * sin(mm*pi*x) .* sin(nn*pi*y);
    fun_phix = @(x,y) mm/(2*pi) * cos(mm*pi*x) .* sin(nn*pi*y);
    fun_phiy = @(x,y) nn/(2*pi) * sin(mm*pi*x) .* cos(nn*pi*y);
elseif srctype == 1
    fun_fext = @(x,y) zeros(size(x));
    fun_phi = @(x,y) 0.1 + 0.2*x + 0.3*y;
    fun_phix = @(x,y) 0.2 * ones(size(x));
    fun_phiy = @(x,y) 0.3 * ones(size(y));
elseif srctype == 2
    fun_fext = @(x,y) -8 .* ones(size(x));
    fun_phi = @(x,y) x.^2 + 2*x.*y + 3 * y.^2;
    fun_phix = @(x,y) 2*x + 2*y;
    fun_phiy = @(x,y) 2*x + 6*y;
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
dilation = 2.1;
% dilation = 3.1;
re = h0 * dilation;

% finite-increment-gradient stabilization
fig_stab_alpha = 0.0;
% fig_stab_alpha = 0.05;
% fig_stab_alpha = 0.5;
% fig_stab_alpha = 0.8;
% fig_stab_alpha = 1.0;

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
    
    pause;
    % return;
end

%
if (1)
    fdpmDriverBuildTri;
    nodeVol = nodeArea;
    pause;
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

if 0
    fdpmDriverPoissonStabHessian;
end

% create sparse K
Ktan0 = SpMatCreateAndClear(Kassem);

disp('End create initial stiffness matrix'); toc;



disp('Begin solve Poisson');

% RHS vector
fext = fun_fext(nodeX,nodeY);
fext = fext .* nodeVol;
% for i = 1:numNodes
    % nneigh = conn(i).numNeigh;
    % ineigh = conn(i).neigh2;
    % ivol = nodeVol(i);
% end

[phi] = fdpmSolve(Ktan0, fext, dirBCDofs,dirBCVals);

disp('End solve Poisson');


if 1
    figure;
    xx = reshape(nodeX, nlen,nlen);
    yy = reshape(nodeY, nlen,nlen);
    ss = reshape(phi, nlen,nlen);
    % surf(xx,yy,ss, 'FaceColor','interp');
    trisurf(tri, nodeX,nodeY,phi, 'FaceColor','interp');
    % zlim([-0.06 0.06]);
    title('solution');
    
    figure;
    zz = fun_phi(xx,yy);
    % surf(xx,yy, zz, 'FaceColor','interp');
    trisurf(tri, nodeX,nodeY, zz(:), 'FaceColor','interp');
    % zlim([-0.06 0.06])
    title('analytical');
end


if 1
    Ktan = Ktan0;
    fdpmDriverShowStiffEigen;
end


if 1
    fdpmDriverShowIntegError;
    vsum0 = v_sum;
end

if 0
    % solve by full integration
    
    gaussord = 7;
    [gaussp, gaussw] = gauss2dTriPoint(gaussord);
    gaussx = gaussp(:,1);
    gaussy = gaussp(:,2);
    gaussn = length(gaussw);
    
    ntri = size(tri, 1);
    
    if 0
        figure;
        plot(nodeX,nodeY,'.');
        axis equal;
        
        hold on;
        
        triplot(tri, nodeX,nodeY);
        
        for k = 1:ntri
            i1 = tri(k,1);
            i2 = tri(k,2);
            i3 = tri(k,3);
            gpx = nodeX(i1) + (nodeX(i2)-nodeX(i1)).*gaussx + (nodeX(i3)-nodeX(i1)).*gaussy;
            gpy = nodeY(i1) + (nodeY(i2)-nodeY(i1)).*gaussx + (nodeY(i3)-nodeY(i1)).*gaussy;
            
            plot(gpx,gpy,'xr');
        end
        
        hold off;
    end
    
    fext2 = zeros(numNodes,1);
    vsum2 = zeros(numNodes,2);
    
    for j = 1:ntri
        i1 = tri(j,1);
        i2 = tri(j,2);
        i3 = tri(j,3);
        x21 = nodeX(i2) - nodeX(i1);
        x31 = nodeX(i3) - nodeX(i1);
        y21 = nodeY(i2) - nodeY(i1);
        y31 = nodeY(i3) - nodeY(i1);
        gpx = nodeX(i1) + x21.*gaussx + x31.*gaussy;
        gpy = nodeY(i1) + y21.*gaussx + y31.*gaussy;
        ga = 0.5 * abs(x21*y31 - y21*x31);
        
        for k = 1:gaussn
            xk = gpx(k);
            yk = gpy(k);
            wk = ga * gaussw(k) / sum(gaussw);
            ok = find(sqrt((xk-nodeX).^2+(yk-nodeY).^2) < re);
            
            [N,Nx,Ny] = mlsRegularShape2D([xk,yk], nodePos(ok,1),nodePos(ok,2), re);
            vsum2(ok,1) = vsum2(ok,1) + wk .* Nx(:);
            vsum2(ok,2) = vsum2(ok,2) + wk .* Ny(:);
            
            dN = [Nx; Ny];
            Ke = dN' * dN .* wk;
            
            SpMatAssemBlockWithDof(Kassem, Ke, ok,ok);
            fext2(ok) = fext2(ok) + fun_fext(xk,yk)*wk .* N(:);
        end
    end
    
    Ktan2 = SpMatCreateAndClear(Kassem);
    [phi2] = fdpmSolve(Ktan2, fext2, dirBCDofs,dirBCVals);
    
    figure; 
    quiver(nodeX(:),nodeY(:),vsum2(:,1),vsum2(:,2));
    axis equal;
    
    figure;
    surf(xx,yy,reshape(phi2,nlen,nlen), 'FaceColor','interp');
    title('full');
end

if 0
    Ktan = Ktan2;
    fdpmDriverShowStiffEigen;
end



return






