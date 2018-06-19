%testPoisson2
% use MLS approx.
% Patch test

clear;



%% geometry
L = 2.0;

%% source function

% srctype = 'sinsin';
srctype = 'linear';
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
    aa = 0.1; bb = 0.2; cc = 0.1;
    fun_fext = @(x,y) zeros(size(x));
    fun_phi = @(x,y) cc + aa*x + bb*y;
    fun_phix = @(x,y) aa * ones(size(x));
    fun_phiy = @(x,y) bb * ones(size(y));
otherwise
    error('unknown srctype: %s', srctype);
end


%


%% fdpm setup

nlen = 2;
% nlen = 4;
h0 = L / nlen;

nlen1 = nlen + 1;


% order
mls_order = 1;
% mls_order = 2;

% influence radius, if use p(2), must > 2
dilation = 1.3;
dilation = 1.6;
% dilation = 2.1;
% dilation = 2.6;
% dilation = 3.1;
% dilation = 4.1;
re = h0 * dilation;


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
if 0
    numNodes = 36;
    nodeX = [ 0.0,0.5,1.0,1.5,2.0,2.0,2.0,2.0,2.0,1.5,1.0,0.5,0.0,0.0,0.0,0.0,0.12,0.15,...
    0.3,0.4,0.55,0.6,0.75,0.8,0.9,1.1,1.2,1.2,1.3,1.4,1.4,1.4,1.5,1.7,1.78,1.9]';
    nodeY = [ 0.0,0.0,0.0,0.0,0.0,0.5,1.0,1.5,2.0,2.0,2.0,2.0,2.0,1.5,1.0,0.5,1.45,0.2,...
    0.98,1.1,1.2,0.8,1.2,1.6,0.4,0.9,0.4,1.9,1.4,0.34,0.9,1.05,0.3,1.8,1.0,0.6]';
    
    numEdges = 16;
    edgePair = [1:16;[2:16,1]].';
end
if 0
    % regular 9 points
    numNodes = 9;
    nodeX = [ 0.0,1.0,2.0,2.0,2.0,1.0,0.0,0.0,1.0 ]';
    nodeY = [ 0.0,0.0,0.0,1.0,2.0,2.0,2.0,1.0,1.0 ]';
    
    numEdges = 8;
    edgePair = [1,2;2,3;3,4;4,5;5,6;6,7;7,8;8,1];
end
if 0
    % regular 8 points and irregular point 9
    numNodes = 9;
    nodeX = [ 0.0,1.0,2.0,2.0,2.0,1.0,0.0,0.0,0.3 ]';
    nodeY = [ 0.0,0.0,0.0,1.0,2.0,2.0,2.0,1.0,0.2 ]';
    
    numEdges = 8;
    edgePair = [1,2;2,3;3,4;4,5;5,6;6,7;7,8;8,1];
end
if 1
    % irregular 15 points
    numNodes = 15;
    nodeX = [ 0.0,0.3,2.0,0.0,1.2,0.8,0.9,1.4,0.5,1.1,2.0,0.0,0.8,1.6,2.0 ]';
    nodeY = [ 0.0,0.0,0.0,0.2,0.4,0.9,1.0,1.2,1.3,1.2,1.8,2.0,2.0,1.9,2.0 ]';
    
    % regular boundary 
    % nodeX = [ 0.0,1.0,2.0,0.0,1.2,0.8,0.9,1.4,0.5,1.1,2.0,0.0,1.0,1.6,2.0 ]';
    % nodeY = [ 0.0,0.0,0.0,1.0,0.4,0.9,1.0,1.2,1.3,1.2,1.0,2.0,2.0,1.9,2.0 ]';
    
    numEdges = 8;
    edgePair = [1,2;2,3;3,11;11,15;15,13;13,12;12,4;4,1];
end

%
nodePos = [nodeX,nodeY];
% 
numDofs = numNodes * 1;


disp('End generate particles & BC');


disp('Begin generate integ points');

numSamples = 0;
sampleX = [];
sampleY = [];
sampleW = [];
if 1
    
    ngauss = 5;
    
    ncell = 16;
    for j = 1:ncell
    for i = 1:ncell
        xc = (i-0.5)*(L/ncell);
        yc = (j-0.5)*(L/ncell);
        ac = (L/ncell)*(L/ncell);
        
        [gn,gp,gw] = gauss2dRectTensorPoint(ngauss,ngauss);
        
        numSamples = numSamples + gn;
        sampleX = [ sampleX; gp(:,1)*(L/ncell)/2+xc ];
        sampleY = [ sampleY; gp(:,2)*(L/ncell)/2+yc ];
        sampleW = [ sampleW; gw(:)/4*ac ];
    end
    end
    
end

sampleR = zeros(numSamples,1);
sampleR(:) = re;
if 1
    for k = 1:numSamples
        kx = sampleX(k);
        ky = sampleY(k);
        
        rx = nodeX - kx;
        ry = nodeY - ky;
        rr = sqrt(rx.^2 + ry.^2);
        [dd,ii] = sort(rr);
        
        sampleR(k) = dd(4) + 2.0e-2;
    end
end


disp('End generate integ points');



if (1)
    % plot fdpm node mesh
    figure;
    plot(nodeX,nodeY,'.', sampleX,sampleY,'+');
    % legend('node');
    % axis([0 L 0 L]);
    axis('equal');
    
    hold on;
    for i = 1:numEdges
        plot(nodeX(edgePair(i,:)),nodeY(edgePair(i,:)),'x-');
    end
    for i = 1:numNodes
        text(nodeX(i),nodeY(i),int2str(i));
    end
    for i = 1:numSamples
        % rectangle('Position',[sampleX(i)-sampleR(i),sampleY(i)-sampleR(i),sampleR(i)*2,sampleR(i)*2],'Curvature',[1,1]);
        % pause
    end
    hold off;
    
    pause;
    % return;
end

%% build triangulation to determine volume
fdpmDriverBuildTri;

% set volume
nodeVol = nodeArea;


%% build MLS
disp('Begin build MLS'); tic;

clear conn1;
for i = 1:numSamples
    
    [neigh,rs,rx,ry,cutoff] = fdpmNeighborhood([sampleX(i),sampleY(i)],nodePos,sampleR(i));
    neigh2 = neigh;
    
    [N,Nx,Ny] = fdpmShapeMLS2(rx(neigh2),ry(neigh2),cutoff, 'mls_order',mls_order, 'mls_volume',nodeVol(neigh2));
    % [N,Nx,Ny] = fdpmShapeMLS3(sampleX(i),sampleY(i),nodeX(neigh2),nodeY(neigh2),cutoff);
    
	% save connection
    conn1(i).numNeigh = numel(neigh2);
	conn1(i).neigh2 = neigh2'; % neighbor list is a row vector
    conn1(i).N = N;
	conn1(i).NX = Nx;
	conn1(i).NY = Ny;
end

disp('End build MLS'); toc;


conn2 = conn1;



%% create initial stiffness matrix
disp('Begin create initial stiffness matrix'); tic;

% assembler
Kassem = SpMatAssem(numDofs);
frhs = zeros(numDofs,1);

for i = 1:numSamples
    nneigh = conn2(i).numNeigh;
    ineigh = conn2(i).neigh2;
    
    ix = sampleX(i);
    iy = sampleY(i);
    ivol = sampleW(i);
    
    N = conn2(i).N;
    dN = [conn2(i).NX; conn2(i).NY];
    Ke = dN' * dN * ivol;
    
    fe = fun_fext(ix,iy);
    
    SpMatAssemBlockWithDof(Kassem, Ke,ineigh,ineigh);
    frhs(ineigh) = frhs(ineigh) + fe*ivol .* N';
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
    
    gnum = 5;
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
            [N] = fdpmShapeMLS2(rx(kneigh),ry(kneigh),cutoff, 'mls_order',mls_order, 'mls_volume',nodeVol(kneigh));
            % [N] = fdpmShapeMLS3(kx,ky,nodeX(kneigh),nodeY(kneigh),cutoff);
            
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

mat = [ Ktan0, Kbnd; Kbnd', sparse(nbndry,nbndry) ];
rhs = [ frhs; fbnd ];
sol = mat \ rhs;
phi = sol(1:numDofs);
lam = sol(numDofs+1:end);


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





