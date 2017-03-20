
func_alpha = @(n,m) sqrt((n+m+1)*(n-m+1));
func_beta = @(n,m) kappa^2 * func_alpha(n,m) / (2*n+1) / (2*n+3);
func_eta = @(n,m) SignNonNeg(m) * sqrt((n-m-1)*(n-m));
func_mu = @(n,m) kappa^2 * func_eta(n,m) / (2*n-1) / (2*n+1);

% Sn,m,l
Smat = zeros(4*nmax1^2);

sind = @(n,l,m) sub2ind(size(Smat), sh_sub2ind(n,m), sh_sub2ind(l,m));

% step 1, n=m=0
for l = 0:2*nmax1-1
	Smat(sind(0,l,0)) = 1/r^l * AdaptModSphBesselK(l,kr) * exp(-kr)/r;
end
% figure; spy(Smat); 
% return

for n = 0:2*nmax1-1
	Smat(sind(n,0,0)) = (-1)^n * Smat(sind(0,n,0));
end
% figure; spy(Smat);
% return

% step 2, m=0
for n = 0:nmax1-2
	for l = n+1:2*nmax1-n-2
		if n == 0
			beta1 = func_beta(l-1,0);
			S1 = Smat(sind(n,l-1,0));
			
			alpha3 = func_alpha(l,0);
			S3 = Smat(sind(n,l+1,0));
			
			alphanew = func_alpha(n,0);
			Snew = -1/alphanew * (beta1*S1 + alpha3*S3);
			Smat(sind(n+1,l,0)) = Snew;
		else
			beta1 = func_beta(l-1,0);
			S1 = Smat(sind(n,l-1,0));
			
			beta2 = func_beta(n-1,0);
			S2 = Smat(sind(n-1,l,0));
			
			alpha3 = func_alpha(l,0);
			S3 = Smat(sind(n,l+1,0));
			
			alphanew = func_alpha(n,0);
			Snew = -1/alphanew * (beta1*S1 + beta2*S2 + alpha3*S3);
			Smat(sind(n+1,l,0)) = Snew;
		end
	end
end
% figure; spy(Smat);
% return


% step 3
for m = 1:nmax1-1
	
	%
	for l = m:2*nmax1-2-m
		
		mu1 = func_mu(l,-m);
		S1 = Smat(sind(m-1,l-1,m-1));
		
		eta2 = func_eta(l+1,m-1);
		S2 = Smat(sind(m-1,l+1,m-1));
		
		etanew = func_eta(m,-m);
		Snew = -1/etanew * (mu1*S1 + eta2*S2);
		
		Smat(sind(m,l,m)) = Snew;
	end
	
	%
	for n = m:nmax1-2
		for l = n+1:2*nmax1-n-3
			if n == m
				beta1 = func_beta(l-1,m);
				S1 = Smat(sind(n,l-1,m));
				
				alpha3 = func_alpha(l,m);
				S3 = Smat(sind(n,l+1,m));
				
				alphanew = func_alpha(n,m);
				Snew = -1/alphanew * (beta1*S1 + alpha3*S3);
				
				Smat(sind(n+1,l,m)) = Snew;
			else
				beta1 = func_beta(l-1,m);
				S1 = Smat(sind(n,l-1,m));
				
				beta2 = func_beta(n-1,m);
				S2 = Smat(sind(n-1,l,m));
				
				alpha3 = func_alpha(l,m);
				S3 = Smat(sind(n,l+1,m));
				
				alphanew = func_alpha(n,m);
				Snew = -1/alphanew * (beta1*S1 + beta2*S2 + alpha3*S3);
				
				Smat(sind(n+1,l,m)) = Snew;
			end
		end
	end
end
% figure; spy(Smat);
% return

% step 4
for n = 1:nmax
	for m = -n:n
		for l = abs(m):nmax1-1
			if n > l
				Smat(sind(n,l,m)) = (-1)^(n+l) * Smat(sind(l,n,m));
			end
		end
	end
end
% figure; spy(Smat);

for n = 1:nmax
	for m = -n:-1
		for l = abs(m):nmax
			Smat(sind(n,l,m)) = Smat(sind(n,l,-m));
		end
	end
end
% figure; spy(Smat);

Smat = Smat(1:npole,1:npole);


