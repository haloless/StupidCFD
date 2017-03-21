
clear;

ImUnit = 1i;

kappa = 1.0;

pa = [ 0; 0; 0 ];
pb = [ 2; 1; 0.5 ];

pab = pa - pb;

[r,theta,phi] = sh_cart2sph(pab(1),pab(2),pab(3));
kr = kappa * r;

epphi = exp(ImUnit*phi);
emphi = exp(-ImUnit*phi);
cost = cos(theta);
sint = sin(theta);

% figure;
% plot3([pa(1),pb(1)],[pa(2),pb(2)],[pa(3),pb(3)],'x-')

nmax = 2;
nmax1 = nmax + 1;
npole = nmax1^2;

func_sign = @(m) sign(m+0.001);
func_anm = @(n,m) sqrt((n+m+1)*(n-m+1)/(2*n+1)/(2*n+3));
func_bnm = @(n,m) func_sign(m) * sqrt((n-m-1)*(n-m)/(2*n-1)/(2*n+1));

Rmat = zeros(4*nmax1^2);

% step 1.
for n = 0:2*nmax1-1
	% m = 0
	irow = sh_sub2ind(n,0);
	
	for s = -n:n
		icol = sh_sub2ind(n,s);
		
		Rmat(irow,icol) = SphHarmY(n,-s,theta,phi);
	end
end

% figure; spy(Rmat);

% step 2.
for m = 0:nmax1-1
	for n = m+2:2*nmax1-m-1
		
		inm = sh_sub2ind(n,m);
		
		% anm = func_anm(n,m);
		bnm = func_bnm(n,m);
		
		for s = -(n-1):n-1
			ins = sh_sub2ind(n,s);
			
			b1 = func_bnm(n,s-1);
			R1 = Rmat(inm,ins-1);
			R1 = 0.5*emphi*(1+cost) * b1 * R1;
			
			b2 = func_bnm(n,-s-1);
			R2 = Rmat(inm,ins+1);
			R2 = 0.5*epphi*(1-cost) * b2 * R2;
			
			% a3 = func_anm(n,s);
			a3 = func_anm(n-1,s);
			R3 = Rmat(inm,ins);
			R3 = sint * a3 * R3;
			
			Rnew = 1/bnm * (R1-R2+R3);
			
			irow = sh_sub2ind(n-1,m+1);
			icol = sh_sub2ind(n-1,s);
			Rmat(irow,icol) = Rnew;
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

% xtmp = rand(npole,1);
% ytmp = zeros(npole,1);
% ztmp = zeros(npole,1);
% for n = 0:nmax
	% for m = -n:n
		% inm = sh_sub2ind(n,m);
		% ynm = 0;
		% for s = -n:n
			% ins = sh_sub2ind(n,s);
			% ynm = ynm + Rmat(inm,ins) * xtmp(ins);
		% end
		% ytmp(inm) = ynm;
	% end
% end
% for n = 0:nmax
	% for m = -n:n
		% inm = sh_sub2ind(n,m);
		% znm = 0;
		% for s = -n:n
			% ins = sh_sub2ind(n,s);
			% znm = znm + conj(Rmat(ins,inm)) * ytmp(ins);
		% end
		% ztmp(inm) = znm;
	% end
% end

%
%

func_alpha = @(n,m) sqrt((n+m+1)*(n-m+1));
func_beta = @(n,m) kappa^2 * func_alpha(n,m) / (2*n+1) / (2*n+3);
func_eta = @(n,m) func_sign(m) * sqrt((n-m-1)*(n-m));
func_mu = @(n,m) kappa^2 * func_eta(n,m) / (2*n-1) / (2*n+1);

% Sn,m,l
Smat = zeros(4*nmax1^2);

% step 1, n=m=0
for l = 0:2*nmax1-1
	irow = sh_sub2ind(0,0);
	icol = sh_sub2ind(l,0);
	Smat(irow,icol) = 1/r^l * AdaptModSphBesselK(l,kr) * exp(-kr)/r;
end
% figure; spy(Smat);

for n = 0:2*nmax1-1
	Smat(sh_sub2ind(n,0),sh_sub2ind(0,0)) = (-1)^n * Smat(sh_sub2ind(0,0),sh_sub2ind(n,0));
end
% figure; spy(Smat);

% step 2, m=0
for n = 0:nmax1-2
	for l = n+1:2*nmax1-n-2
		
		beta1 = func_beta(l-1,0);
		S1 = Smat(sh_sub2ind(n,0),sh_sub2ind(l-1,0));
		
		beta2 = func_beta(n-1,0);
		S2 = Smat(sh_sub2ind(n-1,0),sh_sub2ind(l,0));
		
		alpha3 = func_alpha(l,0);
		S3 = Smat(sh_sub2ind(n,0),sh_sub2ind(l+1,0));
		
		alphanew = func_alpha(n,0);
		Snew = -1/alphanew * (beta1*S1 + beta2*S2 + alpha3*S3);
		Smat(sh_sub2ind(n+1,0),sh_sub2ind(l,0)) = Snew;
	end
end
figure; spy(Smat);
% return
% step 3
for m = 0:nmax1-2
	
	for l = m:2*nmax1-2-m
		
		mu1 = func_mu(l,-m-1);
		S1 = Smat(sh_sub2ind(m,m),sh_sub2ind(l-1,m));
		
		eta2 = func_eta(l+1,m);
		S2 = Smat(sh_sub2ind(m,m),sh_sub2ind(l+1,m));
		
		etanew = func_eta(m+1,-m-1);
		Snew = -1/etanew * (mu1*S1 + eta2*S2);
		
		Smat(sh_sub2ind(m+1,m+1),sh_sub2ind(l,m+1)) = Snew;
	end
	% figure; spy(Smat);
	for n = m:nmax1-2
		for l = n+1:2*nmax1-n-2
			beta1 = func_beta(l-1,m+1);
			S1 = Smat(sh_sub2ind(n,m+1),sh_sub2ind(l-1,m+1));
			
			beta2 = func_beta(n-1,m+1);
			S2 = Smat(sh_sub2ind(n-1,m+1),sh_sub2ind(l,m+1));
			
			alpha3 = func_alpha(l,m+1);
			S3 = Smat(sh_sub2ind(n,m+1),sh_sub2ind(l+1,m+1));
			
			alphanew = func_alpha(n,m+1);
			Snew = -1/alphanew * (beta1*S1 + beta2*S2 + alpha3*S3);
			
			Smat(sh_sub2ind(n+1,m+1),sh_sub2ind(l,m+1)) = Snew;
		end
	end
	% figure; spy(Smat);

end
figure; spy(Smat);
return

% step 4
for n = 0:nmax
	for m = -n:n
		for l = abs(m):nmax
			if n > l
				Smat(sh_sub2ind(n,m),sh_sub2ind(l,m)) = (-1)^(n+l) * Smat(sh_sub2ind(l,m),sh_sub2ind(n,m));
			end
		end
	end
end
figure; spy(Smat);

for n = 0:nmax
	for m = -n:-1
		for l = abs(m):nmax
			Smat(sh_sub2ind(n,m),sh_sub2ind(l,m)) = Smat(sh_sub2ind(n,-m),sh_sub2ind(l,-m));
		end
	end
end
figure; spy(Smat);

Smat = Smat(1:npole,1:npole);



