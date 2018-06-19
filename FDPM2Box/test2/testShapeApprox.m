% Test derivative computation error
%
%
%


clear;


mm = 2; nn = 2;
func_phi = @(x,y) sin(mm*x + nn*y);
func_phix = @(x,y) cos(mm*x + nn*y) * mm;
func_phiy = @(x,y) cos(mm*x + nn*y) * nn;
func_phixx = @(x,y) -sin(mm*x + nn*y) * mm^2;
func_phiyy = @(x,y) -sin(mm*x + nn*y) * nn^2;
func_phixy = @(x,y) -sin(mm*x + nn*y) * mm*nn;


dh = 0.2;
% dh = 0.1;
% dh = 0.05;
% dh = 0.02;
% dh = 0.01;
[Xs,Ys] = ndgrid(-1:dh:1,-1:dh:1);

xs = Xs(:);
ys = Ys(:);

% define a field on current config
phis = func_phi(xs,ys);
phixs = func_phix(xs,ys);
phiys = func_phiy(xs,ys);
phixxs = func_phixx(xs,ys);
phiyys = func_phiyy(xs,ys);
phixys = func_phixy(xs,ys);

if 0
	figure;
	plot(Xs,Ys,'.', xs,ys,'x')
	axis equal
end


numNodes = numel(xs);
nodePos = [xs,ys];

re = dh * 2.1;
% re = dh * 3.1;
fig_stab_alpha = 0.0;

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
    conn(i).wneigh = [w;0];
	conn(i).dNX = dN(1,:);
	conn(i).dNY = dN(2,:);
	conn(i).dNXX = dN(3,:);
	conn(i).dNXY = dN(4,:);
	conn(i).dNYY = dN(5,:);
end
toc;
disp('End build connnectivity');


disp(['dh=',num2str(dh)])
disp(['figstab=',num2str(fig_stab_alpha)])


%% compute dphi/dx
disp('Begin compute g(phi) and h(phi)')
phix = zeros(numNodes,1);
phiy = zeros(numNodes,1);
phixx = zeros(numNodes,1);
phiyy = zeros(numNodes,1);
phixy = zeros(numNodes,1);
for i = 1:numNodes
    ineigh = conn(i).neigh2(:)';
    % ineighpos = nodePos(ineigh,:);
    
    dN = [conn(i).dNX; conn(i).dNY];
    ddN = [conn(i).dNXX; conn(i).dNXY; conn(i).dNYY];
    
    dp = dN * phis(ineigh);
    ddp = ddN * phis(ineigh);
    
    
    phix(i) = dp(1);
    phiy(i) = dp(2);
    phixx(i) = ddp(1);
    phixy(i) = ddp(2);
    phiyy(i) = ddp(3);
end

ephix = fdpmResidError(phix-phixs);
ephiy = fdpmResidError(phiy-phiys);
ephixx = fdpmResidError(phixx-phixxs);
ephixy = fdpmResidError(phixy-phixys);
ephiyy = fdpmResidError(phiyy-phiyys);
[ephix,ephix,ephixx,ephixy,ephiyy]

%%
disp('Begin compute e(phi)');
ephi = zeros(numNodes,1);
for i = 1:numNodes
    ineigh = conn(i).neigh2(:)';
    
    Ng = [conn(i).dNX; conn(i).dNY];
    Nh = [conn(i).dNXX; conn(i).dNXY; conn(i).dNYY];
    
    gi = Ng * phis(ineigh);
    hi = Nh * phis(ineigh);
    
    rx = xs(ineigh) - xs(i);
    ry = ys(ineigh) - ys(i);
    
    ival = [rx,ry]*gi + [0.5*rx.^2,rx.*ry,0.5*ry.^2]*hi + phis(i);
    ierr = (ival-phis(ineigh)).^2 .* conn(i).wneigh;
    ierr = sum(ierr);
    
    ephi(i) = ierr;
end

% fdpmResidError(ephi)
[sum(ephi), mean(ephi)]

%%
if 0
    i = round(numNodes/2);
    ineigh = conn(i).neigh2(:)';
    
    Ng = [conn(i).dNX; conn(i).dNY];
    Nh = [conn(i).dNXX; conn(i).dNXY; conn(i).dNYY];
    
    gi = Ng * phis(ineigh);
    hi = Nh * phis(ineigh);
    
    rx = xs(ineigh) - xs(i);
    ry = ys(ineigh) - ys(i);
    
    ival = [rx,ry]*gi + [0.5*rx.^2,rx.*ry,ry.^2]*hi + phis(i);
    ierr = (ival-phis(ineigh)).^2 .* conn(i).wneigh;
    [sum(ierr),mean(ierr),]
end



return

