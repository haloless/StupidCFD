%testPoisson2
% use MLS approx.
% 

clear;



%% geometry
L = 1.0;

%% source function

srctype = 'sinsin';
% srctype = 'linear';
switch srctype
case 'sinsin'
    mm = 1; nn = 1;
    % mm = 2; nn = 3;
    % f=sin(m*pi*x) * sin(m*pi*y)
    fun_fext = @(x,y) (mm^2+nn^2)/2 * sin(mm*pi*x) .* sin(nn*pi*y);
    fun_phi = @(x,y) 1/(2*pi^2) * sin(mm*pi*x) .* sin(nn*pi*y);
    fun_phix = @(x,y) mm/(2*pi) * cos(mm*pi*x) .* sin(nn*pi*y);
    fun_phiy = @(x,y) nn/(2*pi) * sin(mm*pi*x) .* cos(nn*pi*y);
case 'linear'
    aa = 1; bb = 3; cc = 0;
    fun_fext = @(x,y) zeros(size(x));
    fun_phi = @(x,y) cc + aa*x + bb*y;
    fun_phix = @(x,y) aa * ones(size(x));
    fun_phiy = @(x,y) bb * ones(size(y));
otherwise
    error('unknown srctype: %s', srctype);
end


%


%% fdpm setup

% nlen = 2;
nlen = 10;
% nlen = 20;
% nlen = 40;
h0 = L / nlen;

nlen1 = nlen + 1;


% order
mls_order = 1;
% mls_order = 2;

% influence radius, if use p(2), must > 2
% dilation = 1.8;
dilation = 2.1;
% dilation = 2.5;
% dilation = 3.1;
re = h0 * dilation;

% perturb interal points
perturb = 0.0;
% perturb = 0.4;


% finite-increment-gradient stabilization
fig_stab_alpha = 0.0;

%% generate points and boundary conditions

% generation routine should setup the following data
numNodes = 0;
% particle states
nodeX = [];
nodeY = [];
% BC Dirichlet
dirBCDofs = [];
dirBCVals = [];

disp('Begin generate particles & BC');
if 1
    % generate nodes
    for j = 0:nlen
    for i = 0:nlen
        xx = (i) * h0;
        yy = (j) * h0;
        
        if 1 % perturb internal points, not boundary points
            if i~=0 && i~=nlen && j~=0 && j~=nlen
                rr = perturb;
                xx = xx + (rand(1)-0.5)*h0*rr;
                yy = yy + (rand(1)-0.5)*h0*rr;
            end
        end
        
        numNodes = numNodes + 1;
        
        nn = numNodes;
        nodeX(nn,1) = xx;
        nodeY(nn,1) = yy;
        
        if i==0 || i==nlen || j==0 || j==nlen
            % dirBCDofs(end+1,1) = nn;
            % dirBCVals(end+1,1) = fun_phi(xx,yy);
        end
    end
    end
end

%
nodePos = [nodeX,nodeY];
% 
numDofs = numNodes * 1;

% generate boundary points
numBoundaryNodes = 0;
bndryX = [];
bndryY = [];
bndrySurf = [];
bndryVec = [];

if 1
    ind = reshape(1:numNodes, nlen1,nlen1);
    
    numEdges = 0;
    edgePair = [];
    edgeVec = [];
    for i = 1:nlen
        numEdges = numEdges + 1;
        edgePair(numEdges,:) = [ ind(i,1),ind(i+1,1) ];
        edgeVec(:,numEdges) = [ 0; -1 ];
    end
    for j = 1:nlen
        numEdges = numEdges + 1;
        edgePair(numEdges,:) = [ ind(nlen1,j),ind(nlen1,j+1) ];
        edgeVec(:,numEdges) = [ 1; 0 ];
    end
    for i = nlen:-1:1
        numEdges = numEdges + 1;
        edgePair(numEdges,:) = [ ind(i+1,nlen1),ind(i,nlen1) ];
        edgeVec(:,numEdges) = [ 0; 1 ];
    end
    for j = nlen:-1:1
        numEdges = numEdges + 1;
        edgePair(numEdges,:) = [ ind(1,j+1),ind(1,j) ];
        edgeVec(:,numEdges) = [ -1; 0 ];
    end
    
    for i = 1:numEdges
        i1 = edgePair(i,1);
        i2 = edgePair(i,2);
        x1 = nodeX(i1); y1 = nodeY(i1);
        x2 = nodeX(i2); y2 = nodeY(i2);
        
        ndiv = 1;
        
        for nn = 1:ndiv
            numBoundaryNodes = numBoundaryNodes + 1;
            bndryX(end+1,1) = (nn-0.5) * (x2-x1)/ndiv + x1;
            bndryY(end+1,1) = (nn-0.5) * (y2-y1)/ndiv + y1;
            bndrySurf(end+1,1) = norm([(x2-x1)/ndiv,(y2-y1)/ndiv]);
            bndryVec(:,end+1) = edgeVec(:,i);
        end
    end
end

disp('End generate particles & BC');



if (1)
    % plot fdpm node mesh
    figure;
    plot(nodeX,nodeY,'o', ...
    nodeX(dirBCDofs),nodeY(dirBCDofs),'x');
    legend('node','dir');
    % axis([0 L 0 L]);
    axis('equal');
    
    hold on;
    
    plot(bndryX, bndryY,'.');
    quiver(bndryX,bndryY, bndryVec(1,:)',bndryVec(2,:)');
    hold off;
    
    pause;
    % return;
end

%% build triangulation to determine volume
fdpmDriverBuildTri;

% set volume
nodeVol = nodeArea;


%% build connnectivity

disp('Begin build connnectivity'); tic;
fdpmDriverBuildConn;
disp('End build connnectivity'); toc;

%% build MLS
disp('Begin build MLS'); tic;

clear conn1;
for i = 1:numNodes
    
    [neigh,rs,rx,ry,cutoff] = fdpmNeighborhood(i,nodePos,re);
    neigh2 = [neigh; i]; % particle included in the neighbourhood
    
    [N,Nx,Ny] = fdpmShapeMLS2(rx(neigh2),ry(neigh2),cutoff, 'mls_order',mls_order);
    
	% save connection
    conn1(i).numNeigh = numel(neigh2);
	conn1(i).neigh2 = neigh2'; % neighbor list is a row vector
    conn1(i).N = N;
	conn1(i).NX = Nx;
	conn1(i).NY = Ny;
end

clear connb;
for k = 1:numBoundaryNodes
    [neigh,rs,rx,ry,cutoff] = fdpmNeighborhood([bndryX(k),bndryY(k)], nodePos,re);
    [N,Nx,Ny] = fdpmShapeMLS2(rx(neigh),ry(neigh),cutoff, 'mls_order',mls_order);
    
    connb(k).numNeigh = numel(neigh);
    connb(k).neigh = neigh';
    connb(k).N = N;
    connb(k).NX = Nx;
    connb(k).NY = Ny;
end

disp('End build MLS'); toc;


conn2 = conn1;

% clear conn1;
% for i = 1:numNodes
    
    % [neigh,rs,rx,ry,cutoff] = fdpmNeighborhood(i,nodePos,re);
    % neigh2 = [neigh; i]; % particle included in the neighbourhood
    
    % [N,Nx,Ny,Nxx,Nxy,Nyy] = fdpmShapeMLS2(rx(neigh2),ry(neigh2),cutoff, 'mls_order',2);
    
	% % save connection
    % conn1(i).numNeigh = numel(neigh2);
	% conn1(i).neigh2 = neigh2'; % neighbor list is a row vector
    % conn1(i).N = N;
	% conn1(i).NX = Nx;
	% conn1(i).NY = Ny;
% end


%% build transformation matrix
if 1
    Lmat = sparse(numDofs,numDofs);
    for i = 1:numNodes
        ineigh = conn2(i).neigh2;
        N = conn2(i).N;
        Lmat(ineigh,i) = N.';
    end
    Lmat = Lmat.';
end





%% create initial stiffness matrix
disp('Begin create initial stiffness matrix'); tic;

% assembler
Kassem = SpMatAssem(numDofs);
frhs = zeros(numDofs,1);

fext = fun_fext(nodeX,nodeY);

for i = 1:numNodes
    nneigh = conn2(i).numNeigh;
    ineigh = conn2(i).neigh2;
    ivol = nodeVol(i);
    
    N = conn2(i).N;
    dN = [conn2(i).NX; conn2(i).NY];
    Ke = dN' * dN * ivol;
    
    SpMatAssemBlockWithDof(Kassem, Ke,ineigh,ineigh);
    frhs(ineigh) = frhs(ineigh) + fext(i)*ivol .* N';
    
    % finite increment
    if 1
        eta = 0.0;
        eta = 0.5;
        % eta = 1.0;
        % eta = 2.0;
        hfic = eta * h0;
        
        dH = [conn(i).dNXX; conn(i).dNXY; conn(i).dNXY; conn(i).dNYY];
        Kfic = dH' * dH * ivol * 0.25*hfic^2;
        
        SpMatAssemBlockWithDof(Kassem, Kfic,ineigh,ineigh);
    end
end

% create sparse K
Ktan0 = SpMatCreateAndClear(Kassem);

disp('End create initial stiffness matrix'); toc;



%% boundary part
disp('Begin create boundary matrix'); tic;
if 1
    edgeNode = unique(edgePair(:));
    nbndry = numel(edgeNode);
    Kbnd = sparse(numDofs,nbndry);
    fbnd = zeros(nbndry,1);
    
    gnum = 3;
    [gp,gw] = gauss1dPoint(gnum);
    
    for i = 1:numEdges
        
        ia = edgePair(i,1);
        ib = edgePair(i,2);
        xa = nodeX(ia); ya = nodeY(ia);
        xb = nodeX(ib); yb = nodeY(ib);
        ll = norm([xa,ya]-[xb,yb]);
        
        ind1 = find(edgeNode==ia);
        ind2 = find(edgeNode==ib);
        
        for k = 1:gnum
            kx = (xa+xb)/2 + gp(k)*(xb-xa)/2;
            ky = (ya+yb)/2 + gp(k)*(yb-ya)/2;
            kw = gw(k)*ll/2;
            
            kval = fun_phi(kx,ky);
            
            kn1 = (1-gp(k)) / 2;
            kn2 = (1+gp(k)) / 2;
            
            [kneigh,rs,rx,ry,cutoff] = fdpmNeighborhood([kx,ky], nodePos,re);
            [N] = fdpmShapeMLS2(rx(kneigh),ry(kneigh),cutoff, 'mls_order',mls_order);
            
            Kbnd(kneigh,ind1) = Kbnd(kneigh,ind1) - kw*kn1*N';
            fbnd(ind1) = fbnd(ind1) - kw*kn1*kval;
            
            Kbnd(kneigh,ind2) = Kbnd(kneigh,ind2) - kw*kn2*N';
            fbnd(ind2) = fbnd(ind2) - kw*kn2*kval;
        end
    end
end

disp('End create boundary matrix'); toc;





disp('Begin solve Poisson');

% [phi] = fdpmSolve(Ktan0, frhs, dirBCDofs,dirBCVals);

mat = [ Ktan0, Kbnd; Kbnd', sparse(numBoundaryNodes,numBoundaryNodes) ];
rhs = [ frhs; fbnd ];
sol = mat \ rhs;
phi = sol(1:numDofs);
lam = sol(numDofs+1:end);

phi = Lmat * phi;

disp('End solve Poisson');


if 1
    figure;
    trisurf(tri, nodeX,nodeY,phi, 'FaceColor','interp');
    % zlim([-0.06 0.06]);
    title('solution');
    
    figure;
    zz = fun_phi(nodeX,nodeY);
    trisurf(tri, nodeX,nodeY,zz, 'FaceColor','interp');
    % zlim([-0.06 0.06]);
    title('analytical');
    
    phierr = phi - zz;
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




return



if 0
    % check integration consistency
    gvolx = zeros(numDofs,1);
    gvoly = zeros(numDofs,1);
    gbndx = zeros(numDofs,1);
    gbndy = zeros(numDofs,1);
    
    for i = 1:numNodes
        ineigh = conn2(i).neigh2;
        ivol = nodeVol(i);
        
        Nx = conn2(i).NX';
        Ny = conn2(i).NY';
        
        gvolx(ineigh) = gvolx(ineigh) + ivol.*Nx;
        gvoly(ineigh) = gvoly(ineigh) + ivol.*Ny;
    end
    
    for k = 1:numBoundaryNodes
        kneigh = connb(k).neigh;
        ksurf = bndrySurf(k);
        knvec = bndryVec(:,k);
        
        N = connb(k).N';
        gbndx(kneigh) = gbndx(kneigh) + ksurf*knvec(1).*N;
        gbndy(kneigh) = gbndy(kneigh) + ksurf*knvec(2).*N;
    end
    
    Kassem = SpMatAssem(numDofs);
    for i = 1:numNodes
        ineigh = conn2(i).neigh2;
        ivol = nodeVol(i);
        
        N = conn2(i).N;
        Ki = N' .* ivol;
        
        SpMatAssemBlockWithDof(Kassem, Ki,ineigh,i);
    end
    
    Kcorr = SpMatCreateAndClear(Kassem);
    
    Kcorr = spdiags(nodeVol, 0, numDofs,numDofs) - Kcorr;
    
    if 0
        corrx = Kcorr \ (gbndx-gvolx);
        corry = Kcorr \ (gbndy-gvoly);
    else
        tol = 1.0e-8;
        maxiter = 2000;
        corrx = bicgstab(Kcorr, gbndx-gvolx, tol,maxiter);
        corry = bicgstab(Kcorr, gbndy-gvoly, tol,maxiter);
    end
    
    % apply correction
    for i = 1:numNodes
        ineigh = conn2(i).neigh2;
        ivol = nodeVol(i);
        
        Ni = conn2(i).N;
        Nx = conn2(i).NX;
        Ny = conn2(i).NY;
        
        Nx = Nx - Ni.*corrx(i);
        Ny = Ny - Ni.*corry(i);
        Nx(end) = Nx(end) + corrx(i);
        Ny(end) = Ny(end) + corry(i);
        
        conn2(i).NX = Nx;
        conn2(i).NY = Ny;
    end
    
    % check again
    gvolx = zeros(numDofs,1);
    gvoly = zeros(numDofs,1);
    
    for i = 1:numNodes
        ineigh = conn2(i).neigh2;
        ivol = nodeVol(i);
        
        Nx = conn2(i).NX';
        Ny = conn2(i).NY';
        
        gvolx(ineigh) = gvolx(ineigh) + ivol.*Nx;
        gvoly(ineigh) = gvoly(ineigh) + ivol.*Ny;
    end
    
    gvolx = reshape(gvolx, nlen+1,nlen+1);
    gvoly = reshape(gvoly, nlen+1,nlen+1);
    gbndx = reshape(gbndx, nlen+1,nlen+1);
    gbndy = reshape(gbndy, nlen+1,nlen+1);
    coorx = reshape(corrx, nlen+1,nlen+1);
    coory = reshape(corry, nlen+1,nlen+1);
    
end



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








%%







