% Test derivative computation error
%
%
%


clear;

asmall = 0.25;

func_u = @(x,y) asmall*(x.^2).*y + 0.5*y + 1.5;
func_v = @(x,y) asmall*(y.^3).*(x.^1) + asmall*(x.^2) - 1*x + 2.0;

func_ux = @(x,y) 2*asmall*x.*y;
func_uy = @(x,y) asmall*x.^2 + 0.5;
func_vx = @(x,y) asmall*y.^3 + 2*asmall*x - 1;
func_vy = @(x,y) 3*asmall*x.*(y.^2);

func_uxx = @(x,y) 2*asmall*y;
func_uyy = @(x,y) zeros(size(x));
func_uxy = @(x,y) 2*asmall*x;

func_phi = @(x,y) sin(2*x + 1*y);
func_phix = @(x,y) cos(2*x + 1*y) * 2;
func_phiy = @(x,y) cos(2*x + 1*y) * 1;
func_phixx = @(x,y) -sin(2*x + 1*y) * 4;
func_phiyy = @(x,y) -sin(2*x + 1*y) * 1;
func_phixy = @(x,y) -sin(2*x + 1*y) * 2 * 1;



% dh = 0.2;
dh = 0.1;
% dh = 0.05;
% dh = 0.01;
[Xs,Ys] = ndgrid(-1:dh:1,-1:dh:1);


% current config
xs = Xs(:);
ys = Ys(:);


% define a field on current config
phis = func_phi(xs,ys);
dphidx = func_phix(xs,ys);
dphidy = func_phiy(xs,ys);
dphidxx = func_phixx(xs,ys);
dphidyy = func_phiyy(xs,ys);
dphidxy = func_phixy(xs,ys);

if 0
	figure;
	plot(xs,ys,'.')
	axis equal
    % return
end


numNodes = numel(Xs);
nodePos = [Xs(:),Ys(:)];

re = dh * 1.1;
% re = dh * 2.1;
fig_stab_alpha = 0.0;

if 1
    i = round(numNodes/2);
    hold on;
    plot(xs(i),ys(i),'x');
    hold off;
    
    [neigh,rs,rx,ry,cutoff] = fdpmNeighborhood(nodePos(i,:), nodePos,re);
    
    % [N,Nx,Ny, Nxx,Nxy,Nyy] = fdpmShapeMLS2(rx(neigh),ry(neigh),cutoff);
    [N,Nx,Ny, Nxx,Nxy,Nyy] = fdpmShapeMLS3(xs(i),ys(i),xs(neigh),ys(neigh),cutoff);
    
    dx = 1.0e-5;
    Na = fdpmShapeMLS2(rx(neigh)-dx, ry(neigh), cutoff);
    Nb = fdpmShapeMLS2(rx(neigh)+dx, ry(neigh), cutoff);
    Nxdiff = (Na - Nb) / (dx*2);
    
    % Nxxdiff = (Na - 2*N + Nb) / (dx^2);
    
    % phixxs = func_phixx(xs(i),ys(i));
    % phixx = dot(Nxx, phis(neigh));
    
    return
end

phix = zeros(numNodes,1);
phiy = zeros(numNodes,1);
phixx = zeros(numNodes,1);
phiyy = zeros(numNodes,1);
phixy = zeros(numNodes,1);

for i = 1:numNodes
    [neigh,rs,rx,ry,cutoff] = fdpmNeighborhood(nodePos(i,:), nodePos,re);
    
    [N,Nx,Ny, Nxx,Nxy,Nyy] = fdpmShapeMLS2(rx(neigh),ry(neigh),cutoff);
    
    phix(i) = dot(Nx, phis(neigh));
    phiy(i) = dot(Ny, phis(neigh));
    phixx(i) = dot(Nxx, phis(neigh));
    phiyy(i) = dot(Nyy, phis(neigh));
    phixy(i) = dot(Nxy, phis(neigh));
end

return



disp('Begin build connnectivity');
tic;
clear conn;
for i = 1:numNodes
	
	[neigh,rs,rx,ry,cutoff] = fdpmNeighborhood(i,nodePos,re);
	neigh2 = [neigh; i]; % particle included in the neighbourhood
	
	% weight function
	w = fdpmWeightFunc(rs(neigh),cutoff);
	
	% derivative shape function
	dN = fdpmShapeDer(rx(neigh),ry(neigh),w);
	
	% finite-increment-gradient correction to 1st-derivatives
	if (fig_stab_alpha > 0)
		hX = fig_stab_alpha*dh;
		hY = fig_stab_alpha*dh;
		
		[dN(1,:),dN(2,:)] = fdpmShapeDerFIGStab(dN(1,:),dN(2,:),dN(3,:),dN(4,:),dN(5,:), hX,hY);
	end
	
	% save connection
	conn(i).neigh = neigh;
	conn(i).neigh2 = neigh2;
	conn(i).dNX = dN(1,:);
	conn(i).dNY = dN(2,:);
	conn(i).dNXX = dN(3,:);
	conn(i).dNXY = dN(4,:);
	conn(i).dNYY = dN(5,:);
end
if 0
    % transpose
    for i = 1:numNodes
        ineigh = conn(i).neigh2(:)';
        nneigh = numel(ineigh);
        
        conn(i).dNXt = zeros(1,nneigh);
        conn(i).dNYt = zeros(1,nneigh);
        
        for jj = 1:nneigh
            j = ineigh(jj);
            
            %
            ii = find(conn(j).neigh2 == i);
            
            conn(i).dNXt(jj) = conn(j).dNX(ii);
            conn(i).dNYt(jj) = conn(j).dNY(ii);
        end
    end
end
toc;
disp('End build connnectivity');


uX = zeros(numNodes,1);
uY = zeros(numNodes,1);
vX = zeros(numNodes,1);
vY = zeros(numNodes,1);
F = zeros(2,2,numNodes);
uXX = zeros(numNodes,1);
uYY = zeros(numNodes,1);
uXY = zeros(numNodes,1);

disp('Begin compute F')
for i = 1:numNodes
	ineigh = conn(i).neigh2(:)';
	dNX = conn(i).dNX;
	dNY = conn(i).dNY;
	
	uX(i) = dNX * us(ineigh);
	uY(i) = dNY * us(ineigh);
	vX(i) = dNX * vs(ineigh);
	vY(i) = dNY * vs(ineigh);
	F(:,:,i) = [uX(i), uY(i); vX(i), vY(i)] + eye(2);
    
    uXX(i) = conn(i).dNXX * us(ineigh);
    uYY(i) = conn(i).dNYY * us(ineigh);
    uXY(i) = conn(i).dNXY * us(ineigh);
end

uXerr = fdpmResidError(uX-dudX);
uYerr = fdpmResidError(uY-dudY);
vXerr = fdpmResidError(vX-dvdX);
vYerr = fdpmResidError(vY-dvdY);
disp(['dh=',num2str(dh)])
disp(['figstab=',num2str(fig_stab_alpha)])
[uXerr,uYerr,vXerr,vYerr]

%% compute dphi/dx by dphi/dX
disp('Begin compute grad(phi)')
phix = zeros(numNodes,1);
phiy = zeros(numNodes,1);
for i = 1:numNodes
	ineigh = conn(i).neigh2(:)';
	dNX = conn(i).dNX;
	dNY = conn(i).dNY;
	
	Fi = F(:,:,i);
	
	dphi = Fi' \ [dNX*phis(ineigh); dNY*phis(ineigh)];
	
	phix(i) = dphi(1);
	phiy(i) = dphi(2);
end

phix_err = fdpmResidError(phix-dphidx);
phiy_err = fdpmResidError(phiy-dphidy);
[phix_err,phiy_err]





