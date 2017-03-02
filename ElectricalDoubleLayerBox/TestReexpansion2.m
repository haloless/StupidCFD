
clear all;

ImUnit = 1i;

kappa = 1.0;

pa = [ 0; 0; 0 ];
% pb = [ 2; 1; 0.5 ];
pb = [ 4; 1; 0.5 ];

pab = pa - pb;
% pab = pb - pa;

[r,theta,phi] = sh_cart2sph(pab(1),pab(2),pab(3));
kr = kappa * r;

epphi = exp(ImUnit*phi);
emphi = exp(-ImUnit*phi);
cost = cos(theta);
sint = sin(theta);


ptest = [ 1; 0.5; 0.25 ];
[ra,thetaa,phia] = sh_cart2sph(ptest(1)-pa(1),ptest(2)-pa(2),ptest(3)-pa(3));
[rb,thetab,phib] = sh_cart2sph(ptest(1)-pb(1),ptest(2)-pb(2),ptest(3)-pb(3));



% figure;
% plot3([pa(1),pb(1)],[pa(2),pb(2)],[pa(3),pb(3)],'x-')

% nmax = 2;
nmax = 20;
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


ashy = zeros(npole,1);
bshy = zeros(npole,1);
for n = 0:nmax
	for m = -n:n
		ashy(sh_sub2ind(n,m)) = SphHarmY(n,m,thetaa,phia);
		bshy(sh_sub2ind(n,m)) = SphHarmY(n,m,thetab,phib);
	end
end

aihat = zeros(nmax1,1);
bkhat = zeros(nmax1,1);
for n = 0:nmax
	bkhat(n+1) = exp(-kappa*rb)/rb^(n+1) * AdaptModSphBesselK(n,kappa*rb);
	aihat(n+1) = ra^n * AdaptModSphBesselI(n,kappa*ra);
end

% bcoef = rand(npole,1);
% bcoef = (1.0 ./ (1:npole).^2)';
bcoef = exp(-(1:npole).^2)';
%
Tmat = Rmat' * Smat * Rmat;
% Tmat = Tmat';
lcoef = Tmat * bcoef;


bval = 0;
for n = 0:nmax
	for m = -n:n
		inm = sh_sub2ind(n,m);
		valnm = bkhat(n+1) * bcoef(inm) * bshy(inm);
		bval = bval + valnm;
	end
end

aval = 0;
for n = 0:nmax
	for m = -n:n
		inm = sh_sub2ind(n,m);
		valnm = aihat(n+1) * lcoef(inm) * ashy(inm);
		aval = aval + valnm;
	end
end

