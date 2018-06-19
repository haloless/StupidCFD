
clear;



%% geometry
L = 1.0;

%% source function

% srctype = 0;
srctype = 1;
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
else
    error('Unknown problem');
end


%


%% fdpm setup

nlen = 8;
% nlen = 16;
% nlen = 26;
% nlen = 36;
% nlen = 46;
nlen1 = nlen - 1;
h0 = L / nlen1;


% influence radius, if use p(2), must > 2
% dilation = 1.1;
% dilation = 1.5;
dilation = 1.6;
% dilation = 2.1;
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
        
        if 1 % perturb internal points, not boundary points
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
fdpmDriverBuildConnRegular;
disp('End build connnectivity'); toc;


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
    
    % finite increment
    Kfic = zeros(nneigh,nneigh);
    if 0
        % eta = 0.5;
        eta = 1.0;
        % eta = 2.0;
        % eta = dilation;
        hfic = eta * h0;
        
        dH = [conn(i).dNXX; conn(i).dNXY; conn(i).dNXY; conn(i).dNYY];
        Kfic = dH' * dH * ivol * 0.25*hfic^2;
    end
    
    SpMatAssemBlockWithDof(Kassem, Ke+Kfic,ineigh,ineigh);
end

% create sparse K
Ktan0 = SpMatCreateAndClear(Kassem);

disp('End create initial stiffness matrix'); toc;

if 0
    disp('Begin create stab matrix'); tic;
    
    kstab = 0.0;
    kstab = 0.1;
    % kstab = 0.2;
    % kstab = 0.5;
    % kstab = 1.0;
    % kstab = 2.0;
    if 1
        kstab = kstab / h0^2;
    end
    
    for i = 1:numNodes
        nneigh = conn(i).numNeigh;
        ineigh = conn(i).neigh2;
        wneigh = conn(i).W;
        
        gi = [conn(i).dNX; conn(i).dNY];
        hi = [conn(i).dNXX; conn(i).dNXY; conn(i).dNXY; conn(i).dNYY];
        hi = reshape(hi, 2,2,nneigh);
        % hi(:) = 0;
        
        rx = nodeX(ineigh) - nodeX(i);
        ry = nodeY(ineigh) - nodeY(i);
        mi = [rx, ry, 0.5*rx.^2, rx.*ry, 0.5*ry.^2];
        
        Li = zeros(nneigh);
        if 0
            for jj = 1:nneigh
                j = ineigh(jj);
                rij = [nodeX(j)-nodeX(i); nodeY(j)-nodeY(i)];
                
                for kk = 1:nneigh
                    k = ineigh(kk);
                    Li(jj,kk) = rij'*gi(:,kk) + 0.5*rij'*hi(:,:,kk)*rij;
                    if k == i
                        Li(jj,kk) = Li(jj,kk) + 1;
                    end
                end
            end
        else
            qi = [conn(i).dNX; conn(i).dNY; conn(i).dNXX; conn(i).dNXY; conn(i).dNYY];
            if 1 % use 2nd-order
                Li = mi * qi;
            else % use 1st-order
                Li = mi(:,1:2) * qi(1:2,:);
            end
            Li(:,nneigh) = Li(:,nneigh) + 1;
        end
        
        % part 1
        Ks1 = wneigh * Li - wneigh;
        Ks1 = kstab * Ks1;
        SpMatAssemBlockWithDof(Kassem, Ks1, i,ineigh);
        
        % part 2
        Ks2 = bsxfun(@times, Li, wneigh');
        Ks2 = -kstab * Ks2;
        SpMatAssemBlockWithDof(Kassem, Ks2, ineigh,ineigh);
        
        % part 3
        Ks3 = diag(wneigh);
        Ks3 = kstab * Ks3;
        SpMatAssemBlockWithDof(Kassem, Ks3, ineigh,ineigh);
    end
    
    Kstab = SpMatCreateAndClear(Kassem);
    
    Ktan0 = Ktan0 + Kstab;
    
    disp('End create stab matrix'); toc;
end




disp('Begin solve Poisson');

%
fext = fun_fext(nodeX,nodeY);
fext = fext .* nodeVol;

[phi] = fdpmSolve(Ktan0, fext, dirBCDofs,dirBCVals);

disp('End solve Poisson');


if 1
    figure;
    xx = reshape(nodeX, nlen,nlen);
    yy = reshape(nodeY, nlen,nlen);
    ss = reshape(phi, nlen,nlen);
    surf(xx,yy,ss, 'FaceColor','interp');
    % zlim([-0.06 0.06]);
    title('solution');
    
    figure;
    zz = fun_phi(xx,yy);
    surf(xx,yy, zz, 'FaceColor','interp');
    % zlim([-0.06 0.06])
    title('analytical');
end

if 0
    % check stiffness K matrix eigen structure
    
    figure;
    
    % we copy the K matrix and set Dirichlet BC for boundaries
    Ktan = Ktan0;
    % Ktan(dirBCDofs,dirBCDofs) = 0;
    Ktan(sub2ind([numDofs,numDofs],dirBCDofs,dirBCDofs)) = 1.0e10;
    
    neig = 20;
    % [v,d] = eigs(Ktan, neig);
    [v,d] = eigs(Ktan, neig, 'sm'); % from min eigenvalue
    for ieig = 1:neig
        vv = v(:,ieig);
        surf(xx,yy, reshape(vv,nlen,nlen));
        title(['eig(',int2str(ieig),')=',num2str(d(ieig,ieig))]);
        pause;
    end
end

if 0
    % evaluate energy
    
    % analytical given by variation of 
    % (1/2)*integ(grad_phi.grad_phi dV) - integ(f.phi dV) - integ(grad_phi.n phi dA)
    % note the value of phi vanishes on the boundary so the last boundary integral vanishes
    %
    Eana = integral2(@(x,y) 0.5*(fun_phix(x,y).^2+fun_phiy(x,y).^2) - fun_fext(x,y).*fun_phi(x,y), 0,L,0,L)
    
    Esol = 0;
    for i = 1:numNodes
        nneigh = conn(i).numNeigh;
        ineigh = conn(i).neigh2;
        % wneigh = conn(i).W;
        ivol = nodeVol(i);
        
        dN = [conn(i).dNX; conn(i).dNY];
        gphi = dN * phi(ineigh);
        
        Ei = 0.5*dot(gphi,gphi)*ivol - fext(i)*phi(i);
        Esol = Esol + Ei;
        
        % gphi = 
        % Ei = 
    end
    
    Esol
    
    % Estab = 0.5 * phi' * Kstab * phi
end

if 1
    %
    vsum = zeros(numNodes,2);
    for i = 1:numNodes
        ineigh = conn(i).neigh2;
        ivol = nodeVol(i);
        
        vsum(ineigh,1) = vsum(ineigh,1) + ivol.*conn(i).dNX(:);
        vsum(ineigh,2) = vsum(ineigh,2) + ivol.*conn(i).dNY(:);
    end
    
    figure; 
    quiver(nodeX(:),nodeY(:),vsum(:,1),vsum(:,2));
    axis equal;
end

if 1
    gaussnum = 9;
    [gaussn,gaussp,gaussw] = gauss2dRectTensorPoint(gaussnum, gaussnum);
    gaussx = gaussp(:,1);
    gaussy = gaussp(:,2);
    
    nrect = nlen1;
    % nrect = 4;

    if 0
        figure;
        plot(nodeX,nodeY,'.');
        axis equal;
        
        hold on;
        for j = 1:nrect
        for i = 1:nrect
            xlo = L/nrect * (i-1);
            ylo = L/nrect * (j-1);
            xhi = xlo + L/nrect;
            yhi = ylo + L/nrect;
            
            gpx = (gaussx+1)/2 .*(xhi-xlo) + xlo;
            gpy = (gaussy+1)/2 .*(yhi-ylo) + ylo;
            
            rectangle('Position',[xlo,ylo,xhi-xlo,yhi-ylo]);
            plot(gpx, gpy, 'xr');
        end
        end
        hold off;
    end
    
    
    fext2 = zeros(numNodes,1);
    
    vsum2 = zeros(numNodes,2);
    for j = 1:nrect
    for i = 1:nrect
        xa = L/nrect * (i-1);
        ya = L/nrect * (j-1);
        xb = xa + L/nrect;
        yb = ya + L/nrect;
        
        gpx = (gaussx+1)/2 .*(xb-xa) + xa;
        gpy = (gaussy+1)/2 .*(yb-ya) + ya;
        
        for k = 1:gaussn
            xk = gpx(k);
            yk = gpy(k);
            wk = (xb-xa)*(yb-ya) * gaussw(k) / sum(gaussw);
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

if 1
    vrhs = zeros(numNodes,1);
    for i = 1:numNodes
        nneigh = conn(i).numNeigh;
        ineigh = conn(i).neigh2;
        
        dN = [conn(i).dNX; conn(i).dNY];
        Ke = dN' * dN;
        
        SpMatAssemBlockWithDof(Kassem, Ke,ineigh,ineigh);
        
        for k = 1:nneigh
            j = ineigh(k);
            vrhs(i) = vrhs(i) + dot(vsum2(j,:), dN(:,k));
        end
    end
    vmat = SpMatCreateAndClear(Kassem);
    
end



return






