function [acoef,bcoef] = SolveTwoBodyConstPotent(a1,a2,H,kappa,cutoff, phi1,phi2)


a1k = a1 * kappa;
a2k = a2 * kappa;

% center distance
R = H + a1 + a2;
rk = R * kappa;

% precalc Bessel func.
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

% B
Bmat = zeros(cutoff,cutoff);
for n = 1:cutoff
for m = 1:cutoff
	nn = n - 1;
	mm = m - 1;
	
	Bnm = FuncB(nn,mm,rk);
	
	Bmat(n,m) = Bnm;
end
end

% L,M
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

% identity
Imat = eye(cutoff);

% linear system
mat = [Imat, Lmat; Mmat, Imat];

rhs1 = zeros(cutoff,1);
rhs1(1) = phi1;
rhs2 = zeros(cutoff,1);
rhs2(1) = phi2;
rhs = [ rhs1; rhs2 ];

%
sol = mat \ rhs;

%
acoef = sol(1:cutoff);
bcoef = sol(cutoff+1:end);
acoef = acoef ./ Ka1k;
bcoef = bcoef ./ Ka2k;



return
end

