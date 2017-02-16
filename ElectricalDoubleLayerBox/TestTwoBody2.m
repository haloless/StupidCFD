
clear all;

a1 = 1.0;
a2 = 1.0;

% 
% kappa = 5.0;
% kappa = 4.0;
% kappa = 1.0;
% kappa = 0.5;
kappa = 0.1;
% 
a1k = a1 * kappa;
a2k = a2 * kappa;


% separation
% H = 2.0;
H = 1.0;
% H = 0.5;
% H = 0.2;
% center distance
R = H + a1 + a2;

phi1 = 1.0;
phi2 = 1.0;


cutoff = 20;
cutoff2 = cutoff;


Ka1k = zeros(cutoff,1);
Ia1k = zeros(cutoff,1);
Ka2k = zeros(cutoff,1);
Ia2k = zeros(cutoff,1);
for n = 1:cutoff
	nn = n - 1;
	Ka1k(n) = ModSphBesselK(nn, a1k);
	Ia1k(n) = ModSphBesselI(nn, a1k);
	Ka2k(n) = ModSphBesselK(nn, a2k);
	Ia2k(n) = ModSphBesselI(nn, a2k);
end

Bmat = zeros(cutoff,cutoff);
for n = 1:cutoff
for m = 1:cutoff
	nn = n - 1;
	mm = m - 1;
	
	Bnm = FuncB(nn,mm,kappa*R);
	
	Bmat(n,m) = Bnm;
end
end

Imat = eye(cutoff);
Lmat = zeros(cutoff,cutoff);
Mmat = zeros(cutoff,cutoff);
for j = 1:cutoff
for n = 1:cutoff
	jj = j - 1;
	nn = n - 1;
	
	Bnj = Bmat(n,j);
	Ljn = (2*jj+1) * Bnj * Ia1k(j) / Ka2k(n);
	Mjn = (2*jj+1) * Bnj * Ia2k(j) / Ka1k(n);
	
	Lmat(j,n) = Ljn;
	Mmat(j,n) = Mjn;
end
end

mat = [Imat, Lmat; Mmat, Imat];

rhs1 = zeros(cutoff,1);
rhs1(1) = phi1;
rhs2 = zeros(cutoff,1);
rhs2(1) = phi2;
rhs = [ rhs1; rhs2 ];

sol = mat \ rhs;

acoef = sol(1:cutoff);
bcoef = sol(cutoff+1:end);
acoef = acoef ./ Ka1k;
bcoef = bcoef ./ Ka2k;


[xs,ys] = ndgrid(-a1*2:0.1:R+a2*2,0:0.1:a1*2);
nx = size(xs,1);
ny = size(xs,2);
r1s = sqrt((xs-0).^2 + ys.^2);
r2s = sqrt((xs-R).^2 + ys.^2);
t1s = atan2(ys,xs-0);
t2s = atan2(ys,R-xs);
ct1s = cos(t1s);
ct2s = cos(t2s);

ok = ones(nx,ny);
ok(r1s<a1) = 0;
ok(r2s<a2) = 0;

phis = zeros(nx,ny);

for i = 1:nx
for j = 1:ny
	if ok == 0 
		continue;
	end
	
	xx = xs(i,j);
	yy = ys(i,j);
	r1 = r1s(i,j);
	r2 = r2s(i,j);
	r1k = r1 * kappa;
	r2k = r2 * kappa;
	t1 = t1s(i,j);
	t2 = t2s(i,j);
	ct1 = ct1s(i,j);
	ct2 = ct2s(i,j);
	
	phi = 0;
	for n = 1:cutoff2
		nn = n - 1;
		
		an = acoef(n);
		bn = bcoef(n);
		
		kn1 = ModSphBesselK(nn,r1k);
		kn2 = ModSphBesselK(nn,r2k);
		pn1 = LegendrePoly(nn,ct1);
		pn2 = LegendrePoly(nn,ct2);
		
		un = an*kn1*pn1 + bn*kn2*pn2;
		
		phi = phi + un;
	end
	
	phis(i,j) = phi;
end
end

if 1
phis(ok==0) = nan;
figure; contourf(xs,ys,phis); axis equal
end

