
epphi = exp(ImUnit*phi);
emphi = exp(-ImUnit*phi);
cost = cos(theta);
sint = sin(theta);

func_anm = @(n,m) sqrt((n+m+1)*(n-m+1)/(2*n+1)/(2*n+3));
func_bnm = @(n,m) SignNonNeg(m) * sqrt((n-m-1)*(n-m)/(2*n-1)/(2*n+1));

Rmat = zeros(4*nmax1^2);

rind = @(n,m,s) sub2ind(size(Rmat), sh_sub2ind(n,m), sh_sub2ind(n,s));


% step 1., m=0
for n = 0:2*nmax1-1
	for s = -n:n
		Rmat(rind(n,0,s)) = SphHarmY(n,-s,theta,phi);
	end
end

% figure; spy(Rmat);

% step 2.
for m = 0:nmax1-1
	for n = m+2:2*nmax1-m-1
		
		% anm = func_anm(n,m);
		bnm = func_bnm(n,m);
		
		for s = -n+1:n-1
			
			b1 = func_bnm(n,s-1);
			R1 = Rmat(rind(n,m,s-1));
			R1 = 0.5*emphi*(1+cost) * b1 * R1;
			
			b2 = func_bnm(n,-s-1);
			R2 = Rmat(rind(n,m,s+1));
			R2 = 0.5*epphi*(1-cost) * b2 * R2;
			
			a3 = func_anm(n-1,s);
			R3 = Rmat(rind(n,m,s));
			R3 = sint * a3 * R3;
			
			% Rnew = 1/bnm * (R1-R2+R3);
			% Rnew = -exp(ImUnit*pi)/bnm * (R1-R2+R3);
			Rnew = -exp(ImUnit*angchi)/bnm * (R1-R2+R3);
			
			Rmat(rind(n-1,m+1,s)) = Rnew;
		end
	end
	% figure; spy(Rmat);
end

for n = 0:nmax
	for m = -n:-1
		for s = -n:n
			Rmat(sh_sub2ind(n,m),sh_sub2ind(n,s)) = conj(Rmat(sh_sub2ind(n,-m),sh_sub2ind(n,-s)));
		end
	end
	% figure; spy(Rmat);
end

Rmat = Rmat(1:npole,1:npole);


