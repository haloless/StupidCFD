
clear all;

addpath('../../SphereGridBox');

mu = 1.0 / (6*pi);
g8 = 8*pi*mu;
g6 = 6*pi*mu;

asph = 1;
xsph = [ 0.0, 0.0, 0.0 ]';

dh = asph / 4;
np = SpiralPointsEstimNum(dh/asph)

pos = SpiralPoints(asph,xsph,np);

xi = 2.5 / dh;
alpha = 0;

rcut = 1000.0 * dh;
% rcut = 2.0 * dh;

% figure;
% plot3(pos(1,:),pos(2,:),pos(3,:),'x');

Mmat = zeros(np*3,np*3);
for i = 1:np
for j = 1:np
	i3 = (i-1) * 3;
	j3 = (j-1) * 3;
	
	xij = pos(:,i) - pos(:,j);
	rij = sqrt(xij' * xij);
	if i == j
		rij = 0;
	end
	
	if rij > rcut
		continue;
	end
	
	gij = FuncGReg(rij,xij,alpha,xi);
	gij = gij ./ g8;
	
	Mmat(i3+1:i3+3,j3+1:j3+3) = gij;
end
end

Rmat = inv(Mmat);

Smat = zeros(6,np*3);
for i = 1:np
	i3 = (i-1)*3;
	xx = pos(:,i) - xsph;
	Smat(1:3,i3+1:i3+3) = eye(3);
	Smat(4:6,i3+1:i3+3) = [ 0,-xx(3),xx(2); xx(3),0,-xx(1); -xx(2),xx(1),0];
end

Rcomp = Smat * Rmat * Smat' ./ (g6*asph)







