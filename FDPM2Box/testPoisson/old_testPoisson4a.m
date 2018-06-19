%testPoisson4a
% use a pseudo gauss quadrature
% gauss points are associated with nodes
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

%%
gaussnum = 1;
% gaussnum = 2;
% gaussnum = 3;
[gausspos,gausswgt] = gauss1dPoint(gaussnum);

% 2d gauss
gausstype = gaussnum^2;
[gaussxx,gaussyy] = ndgrid(gausspos,gausspos);
gaussww = gausswgt * gausswgt';
% 
gaussxx = gaussxx(:);
gaussyy = gaussyy(:);
gaussww = gaussww(:);



%% fdpm setup

% nlen = 2;
% nlen = 8;
nlen = 10;
% nlen = 16;
% nlen = 20;
% nlen = 26;
% nlen = 36;
% nlen = 46;
h0 = L / nlen;

nlen1 = nlen;
% nlen1 = nlen + 1;

% influence radius, if use p(2), must > 2
dilation = 2.1;
% dilation = 3.1;
re = h0 * dilation;

% finite-increment-gradient stabilization
fig_stab_alpha = 0.0;

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
    % generate nodes
    for j = 1:nlen
    for i = 1:nlen
        xx = (i-0.5) * h0;
        yy = (j-0.5) * h0;
        
        vv = h0^2;
        
        if 0 % perturb internal points, not boundary points
            if i~=1 && i~=nlen && j~=1 && j~=nlen
                rr = 0.1;
                xx = xx + (rand(1)-0.5)*h0*rr;
                yy = yy + (rand(1)-0.5)*h0*rr;
            end
        end
        
        numNodes = numNodes + 1;
        
        nn = numNodes;
        nodeX(nn,1) = xx;
        nodeY(nn,1) = yy;
        nodeVol(nn,1) = vv;
        
        if i==1 || i==nlen || j==1 || j==nlen
            % dirBCDofs(end+1,1) = nn;
            % dirBCVals(end+1,1) = fun_phi(xx,yy);
        end
    end
    end
end
if 0
    % generate nodes
    for j = 0:nlen
    for i = 0:nlen
        xx = (i) * h0;
        yy = (j) * h0;
        
        vv = h0^2;
        % if i==1 || i==nlen; vv = vv / 2; end;
        % if j==1 || j==nlen; vv = vv / 2; end;
        
        if 0 % perturb internal points, not boundary points
            if i~=0 && i~=nlen && j~=0 && j~=nlen
                rr = 0.4;
                xx = xx + (rand(1)-0.5)*h0*rr;
                yy = yy + (rand(1)-0.5)*h0*rr;
            end
        end
        
        numNodes = numNodes + 1;
        
        nn = numNodes;
        nodeX(nn,1) = xx;
        nodeY(nn,1) = yy;
        nodeVol(nn,1) = vv;
        
        if i==0 || i==nlen || j==0 || j==nlen
            dirBCDofs(end+1,1) = nn;
            dirBCVals(end+1,1) = fun_phi(xx,yy);
        end
    end
    end
end

%
numBnds = 0;
bndX = [];
bndY = [];

% place NN control points on each side
nn = nlen*1;
% nn = nlen*2;
% nn = nlen*4;
% nn = nlen*8;

xx = ((1:nn).'-0.5) * (L/nn);
yy = 0;
numBnds = numBnds + nn;
bndX(end+1:end+nn,1) = xx;
bndY(end+1:end+nn,1) = yy;

xx = ((1:nn).'-0.5) * (L/nn);
yy = L;
numBnds = numBnds + nn;
bndX(end+1:end+nn,1) = xx;
bndY(end+1:end+nn,1) = yy;

xx = 0;
yy = ((1:nn).'-0.5) * (L/nn);
numBnds = numBnds + nn;
bndX(end+1:end+nn,1) = xx;
bndY(end+1:end+nn,1) = yy;

xx = L;
yy = ((1:nn).'-0.5) * (L/nn);
numBnds = numBnds + nn;
bndX(end+1:end+nn,1) = xx;
bndY(end+1:end+nn,1) = yy;


disp('End generate particles & BC');

% 
nodePos = [nodeX,nodeY];

% 
numDofs = numNodes * 1;

%% generate quadrature points

disp('Begin generate integ points');
if 1
    numPatch = numNodes;
    numInteg = numNodes * gausstype;
    integX = zeros(gausstype,numNodes);
    integY = zeros(gausstype,numNodes);
    integW = zeros(gausstype,numNodes);
    for i = 1:numNodes
        ix = nodeX(i);
        iy = nodeY(i);
        
        xx = (h0/2).*gaussxx + ix;
        yy = (h0/2).*gaussyy + iy;
        
        ww = (h0^2)/sum(gaussww) .* gaussww;
        
        integX(:,i) = xx(:);
        integY(:,i) = yy(:);
        integW(:,i) = ww(:);
    end
end
if 0
    numPatch = nlen * nlen;
    numInteg = numPatch * gausstype;
    integX = zeros(gausstype,numPatch);
    integY = zeros(gausstype,numPatch);
    integW = zeros(gausstype,numPatch);
    for j = 1:nlen
    for i = 1:nlen
        ix = (i-0.5)*h0;
        iy = (j-0.5)*h0;
        
        [xx,yy] = ndgrid(gausspos,gausspos);
        xx = (h0/2).*xx + ix;
        yy = (h0/2).*yy + iy;
        
        ww = gausswgt * gausswgt';
        ww = (h0^2)/sum(ww(:)) .* ww;
        
        idx = i + (j-1)*nlen;
        integX(:,idx) = xx(:);
        integY(:,idx) = yy(:);
        integW(:,idx) = ww(:);
    end
    end
end
disp('End generate integ points');


if (1)
    % plot fdpm node mesh
    figure;
    plot(nodeX,nodeY,'o', ...
    integX(:),integY(:),'.',...
    bndX(:),bndY(:),'+');
    legend('node','integ','bnd');
    axis([0-0.2 L+0.2 0-0.2 L+0.2]);
    % axis([0 L 0 L]);
    axis('equal');
    
    hold on;
    rectangle('Position',[0,0,L,L]);
    
    hold off;
    
    pause;
    % return;
end


%% build connnectivity

disp('Begin build connnectivity'); tic;
fdpmDriverBuildConn;
disp('End build connnectivity'); toc;


%% build MLS

disp('Begin build integ MLS'); tic;
clear conni;
clear connp;


nodeOrder = zeros(numNodes,1);
nodeOrder(:) = 2;

bndOrder = zeros(numBnds,1);
bndOrder(:) = 1;

if 1
     for k = 1:numBnds
        xk = bndX(k); yk = bndY(k);
        cutk = re;
        
        %
        [neigh,rs,rx,ry,cutoff] = fdpmNeighborhood([xk,yk],nodePos,cutk);
        
        nodeOrder(neigh) = 1;
    end
end


if 0
    for i = 1:numPatch
        for l = 1:gausstype
            
            xx = integX(l,i); yy = integY(l,i);
            [neigh,rs,rx,ry,cutoff] = fdpmNeighborhood([xx,yy],nodePos,re);
            
            % [N,Nx,Ny] = fdpmShapeMLS2(rx(neigh),ry(neigh),cutoff);
            [N,Nx,Ny] = fdpmShapeMLS3(xx,yy,nodeX(neigh),nodeY(neigh),cutoff);
            
            % save connection
            conni(l,i).numNeigh = numel(neigh);
            conni(l,i).neigh = neigh'; % neighbor list is a row vector
            conni(l,i).N = N;
            conni(l,i).NX = Nx;
            conni(l,i).NY = Ny;
        end
    end
end

if 1
    for i = 1:numNodes
        %
        ix = nodeX(i); iy = nodeY(i);
        [neigh,rs,rx,ry,cutoff] = fdpmNeighborhood([ix,iy],nodePos,re);
        
        [N,Nx,Ny, ~,~,~, ah] = fdpmShapeMLS2(rx(neigh),ry(neigh),cutoff, 'mls_order',nodeOrder(i));
        % [N,Nx,Ny] = fdpmShapeMLS3(xx,yy,nodeX(neigh),nodeY(neigh),cutoff);
        
        if 0
            connp(i).numNeigh = numel(neigh);
            connp(i).neigh = neigh';
            connp(i).N = N;
            connp(i).Nx = Nx;
            connp(i).Ny = Ny;
        end
        
        for l = 1:gausstype
            
            xx = integX(l,i) - ix;
            yy = integY(l,i) - iy;
            
            ph = [ 1, xx, yy, 1/2*xx^2, xx*yy, 1/2*yy^2 ];
            phx = [ 0, 1, 0, xx, yy, 0 ];
            phy = [ 0, 0, 1, 0, xx, yy ];
            
            mh = size(ah,1);
            ph = ph(1:mh); phx = phx(1:mh); phy = phy(1:mh);
            
            N = ph * ah;
            Nx = phx * ah;
            Ny = phy * ah;
            
            conni(l,i).numNeigh = numel(neigh);
            conni(l,i).neigh = neigh'; % neighbor list is a row vector
            conni(l,i).N = N;
            conni(l,i).NX = Nx;
            conni(l,i).NY = Ny;
        end
    end
end
disp('End build integ MLS'); toc;



%% integration consistency

disp('Begin check integration consistency');
if 0
    % check integration consistency
    
    % volume part
    gvolx = zeros(numDofs,1);
    gvoly = zeros(numDofs,1);
    
    for i = 1:numNodes
        ineigh = conn1(i).neigh2;
        ivol = nodeVol(i);
        
        Nx = conn1(i).dNX';
        Ny = conn1(i).dNY';
        
        gvolx(ineigh) = gvolx(ineigh) + ivol.*Nx;
        gvoly(ineigh) = gvoly(ineigh) + ivol.*Ny;
    end
    
    % boundary part
    gbndx = zeros(numDofs,1);
    gbndy = zeros(numDofs,1);
    
    for jj = 1:nlen
    for ii = 1:nlen
        k = ii + (jj-1)*nlen;
        
        dbx = 0;
        dby = 0;
        if ii==1 || ii==nlen
            dby = h0;
            if jj==1 || jj==nlen
                % dby = h0/2;
            end
        end
        if jj==1 || jj==nlen
            dbx = h0;
            if ii==1 || ii==nlen
                % dbx = h0/2;
            end
        end
        
        gbx = 0;
        gby = 0;
        if ii == 1
            gbx = gbx - dby;
        elseif ii == nlen
            gbx = gbx + dby;
        end
        if jj == 1
            gby = gby - dbx;
        elseif jj == nlen
            gby = gby + dbx;
        end
        
        gbndx(k) = gbndx(k) + gbx;
        gbndy(k) = gbndy(k) + gby;
    end
    end
    
    % correction matrix
    Kassem = SpMatAssem(numDofs);
    for i = 1:numNodes
        nneigh = conn1(i).numNeigh;
        ineigh = conn1(i).neigh2;
        ivol = nodeVol(i);
        
        hi = [conn1(i).dNXX; conn1(i).dNXY; conn1(i).dNXY; conn1(i).dNYY];
        
        for nn = 1:nneigh
            j = ineigh(nn);
            
            % hij = hi(:,:,nn);
            hij = sum(hi(:,nn));
            
            SpMatAssemBlockWithDof(Kassem, hij*ivol, j,i);
        end
    end
    
    Kcorr = SpMatCreateAndClear(Kassem);
    
    rcorr = [gbndx-gvolx, gbndy-gvoly];
    
    if 0
        gcorr = Kcorr \ rcorr;
    else
        tol = 1.0e-5;
        maxiter = 2000;
        corrx = bicgstab(Kcorr, rcorr(:,1), tol,maxiter);
        corry = bicgstab(Kcorr, rcorr(:,2), tol,maxiter);
        
        figure;
        plot(rv1/norm(rcorr));
    end
    
    
    % % apply correction
    % for i = 1:numNodes
        % ineigh = conn2(i).neigh2;
        % ivol = nodeVol(i);
        
        % Ni = conn2(i).N;
        % Nx = conn2(i).NX;
        % Ny = conn2(i).NY;
        
        % Nx = Nx - Ni.*corrx(i);
        % Ny = Ny - Ni.*corry(i);
        % Nx(end) = Nx(end) + corrx(i);
        % Ny(end) = Ny(end) + corry(i);
        
        % conn2(i).NX = Nx;
        % conn2(i).NY = Ny;
    % end
    
    % % check again
    % gvolx = zeros(numDofs,1);
    % gvoly = zeros(numDofs,1);
    
    % for i = 1:numNodes
        % ineigh = conn2(i).neigh2;
        % ivol = nodeVol(i);
        
        % Nx = conn2(i).NX';
        % Ny = conn2(i).NY';
        
        % gvolx(ineigh) = gvolx(ineigh) + ivol.*Nx;
        % gvoly(ineigh) = gvoly(ineigh) + ivol.*Ny;
    % end
    
    gvolx = reshape(gvolx, nlen,nlen);
    gvoly = reshape(gvoly, nlen,nlen);
    gbndx = reshape(gbndx, nlen,nlen);
    gbndy = reshape(gbndy, nlen,nlen);
    coorx = reshape(corrx, nlen,nlen);
    coory = reshape(corry, nlen,nlen);
    
    return
end
disp('End check integration consistency');





%% create initial stiffness matrix
disp('Begin create initial stiffness matrix'); tic;

% assembler
Kassem = SpMatAssem(numDofs);

% fext = fun_fext(nodeX,nodeY);
% frhs = fext .* nodeVol;
frhs = zeros(numDofs,1);

for i = 1:numPatch
    
    for l = 1:gausstype
        nneigh = conni(l,i).numNeigh;
        ineigh = conni(l,i).neigh;
        ivol = integW(l,i);
        ix = integX(l,i);
        iy = integY(l,i);
        
        N = conni(l,i).N;
        dN = [conni(l,i).NX; conni(l,i).NY];
        
        %
        Ke = dN' * dN * ivol;
        
        SpMatAssemBlockWithDof(Kassem, Ke, ineigh,ineigh);
        
        %
        frhs(ineigh) = frhs(ineigh) + ivol * fun_fext(ix,iy) * N.';
    end
    
    % finite increment stab.
    if 0
        Kfic = zeros(nneigh,nneigh);
        eta = 0.0;
        % eta = 0.5;
        % eta = 1.0;
        % eta = 2.0;
        % eta = dilation;
        hfic = eta * h0;
        
        dH = [conn1(i).dNXX; conn1(i).dNXY; conn1(i).dNXY; conn1(i).dNYY];
        Kfic = dH' * dH * ivol * 0.25*hfic^2;
        
        % assemble
        SpMatAssemBlockWithDof(Kassem, Kfic,ineigh,ineigh);
    end
    
    
end

% create sparse K
Ktan0 = SpMatCreateAndClear(Kassem);

disp('End create initial stiffness matrix'); toc;


if 1
    kpenal = 1.0e6;
    % kpenal = 1.0e8;
    kpenal = kpenal * (L*4/numBnds);
    
    Kassem = SpMatAssem(numDofs);
    fpenal = zeros(numDofs,1);
    
    for k = 1:numBnds
        xk = bndX(k); yk = bndY(k);
        
        % mls = 1; cutk = re;
        % mls = 2; cutk = re + h0;
        
        %
        [neigh,rs,rx,ry,cutoff] = fdpmNeighborhood([xk,yk],nodePos,re);
        
        [N,Nx,Ny] = fdpmShapeMLS2(rx(neigh),ry(neigh),cutoff, 'mls_order',bndOrder(k));
        % [N,Nx,Ny] = fdpmShapeMLS3(xx,yy,nodeX(neigh),nodeY(neigh),cutoff);
        
        Ke = kpenal * (N'*N);
        SpMatAssemBlockWithDof(Kassem, Ke, neigh,neigh);
        
        phik = fun_phi(xk,yk);
        fpenal(neigh) = fpenal(neigh) + kpenal*phik*N';
    end
    
    Kpenal = SpMatCreateAndClear(Kassem);
    
    Ktan0 = Ktan0 + Kpenal;
    frhs = frhs + fpenal;
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




disp('Begin solve Poisson');


[phi] = fdpmSolve(Ktan0, frhs, dirBCDofs,dirBCVals);

disp('End solve Poisson');


if 1
    xx = reshape(nodeX,nlen1,nlen1);
    yy = reshape(nodeY,nlen1,nlen1);
    
    figure;
    pp = reshape(phi,nlen1,nlen1);
    % tri = delaunay(nodeX,nodeY);
    % trisurf(tri, nodeX,nodeY,phi, 'FaceColor','interp');
    surf(xx,yy,pp, 'FaceColor','interp');
    % zlim([-0.06 0.06]);
    title('solution');
    
    figure;
    zz = fun_phi(xx,yy);
    % trisurf(tri, nodeX,nodeY,zz, 'FaceColor','interp');
    % zlim([-0.06 0.06]);
    surf(xx,yy,zz, 'FaceColor','interp');
    title('analytical');
    
    phierr = phi(:) - zz(:);
end





if 0
    % check stiffness K matrix eigen structure
    
    figure;
    
    % we copy the K matrix and set Dirichlet BC for boundaries
    Ktan = Ktan0;
    % Ktan(sub2ind([numDofs,numDofs],dirBCDofs,dirBCDofs)) = 1.0e10;
    
    neig = 20;
    % [v,d] = eigs(Ktan, neig);
    [v,d] = eigs(Ktan, neig, 'sm'); % from min eigenvalue
    for ieig = 1:neig
        vv = v(:,ieig);
        vv = real(vv);
        surf(xx,yy, reshape(vv,nlen1,nlen1));
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
    
    Estab = 0.5 * phi' * Kstab * phi
end



return






