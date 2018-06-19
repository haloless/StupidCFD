
clear;

kappa = 2.0;

pa = [ 0, 0, 0 ]';
% pb = [ 4, 0, 0 ]';
% pb = [ 2.5, 0, 0 ]';
% pb = [ 0, 0, 2.5 ]';
pb = [ 0, 2.5, 0 ]';
% pb = [ 3, 0, 0 ]';
pb = rand(3,1);
% pb = pb./norm(pb) * 3.0;
pb = pb./norm(pb) * 2.1;

% pb = [ 40, 0, 0 ]';

ra = 1.0;
rb = 1.0;

nmax = 6;
nmax1 = nmax + 1;
npole = nmax1^2;
disp(['nmax=',int2str(nmax), '; npole=',int2str(npole)]);

func_reg = @(n,r) r^n * AdaptModSphBesselI(n,kappa*r);
func_sing = @(n,r) exp(-kappa*r) / r^(n+1) * AdaptModSphBesselK(n,kappa*r);



Tab = ReExpMat(pb,pa, kappa,nmax);
Tba = ReExpMat(pa,pb, kappa,nmax);

khata = zeros(npole,1);
ihata = zeros(npole,1);
khatb = zeros(npole,1);
ihatb = zeros(npole,1);
for n = 0:nmax
	kna = func_sing(n,ra);
	knb = func_sing(n,rb);
	ina = func_reg(n,ra);
	inb = func_reg(n,rb);
	for m = -n:n
		inm = sh_sub2ind(n,m);
		khata(inm) = kna;
		ihata(inm) = ina;
		khatb(inm) = knb;
		ihatb(inm) = inb;
	end
end

La = diag(khata);
Ma = diag(ihata) * conj(Tab);
Lb = diag(khatb);
Mb = diag(ihatb) * conj(Tba);

rhsa = zeros(npole,1);
rhsb = zeros(npole,1);
rhsa(1) = 1.0;
rhsb(1) = 1.0;

A = [La,Ma; Mb,Lb];
rhs = [rhsa;rhsb];
sol = A \ rhs;

coefa = sol(1:npole);
coefb = sol(npole+1:end);


if 0
	% create plot
	
	[xg,yg] = ndgrid(-2:0.2:6,0:0.2:4);
	ug = zeros(size(xg));
	
	for ig = 1:size(ug,1)
	for jg = 1:size(ug,2)
		ptest = [xg(ig,jg),yg(ig,jg),0]';
		
		[rada,thetaa,phia] = sh_cart2sph(ptest-pa);
		[radb,thetab,phib] = sh_cart2sph(ptest-pb);
		
		if rada>=ra && radb>=rb
			shya = zeros(npole,1);
			shyb = zeros(npole,1);
			for n = 0:nmax
			for m = -n:n
				inm = sh_sub2ind(n,m);
				shya(inm) = func_sing(n,rada) * SphHarmY(n,m,thetaa,phia);
				shyb(inm) = func_sing(n,radb) * SphHarmY(n,m,thetab,phib);
			end
			end
			vala = sum(coefa.*shya);
			valb = sum(coefb.*shyb);
			val = vala + valb;
			
			ug(ig,jg) = val;
		else
			ug(ig,jg) = nan;
		end
	end
	end
	
	figure;
	contourf(xg,yg,ug);
	axis equal;
	colorbar;
end

if 1
	% interaction energy
	[khats,kprimes] = ReExpFuncSingular(nmax,ra,kappa);
	[ihats,iprimes] = ReExpFuncRegular(nmax,ra,kappa);
	
	as = coefa;
	bs = conj(Tab) * coefb;
	
	ene = 0;
	for n = 0:nmax
		for m = -n:n
			inm = sh_sub2ind(n,m);
			c1 = as(inm)*khats(n+1) + bs(inm)*ihats(n+1);
			
			inm = sh_sub2ind(n,-m);
			c2 = as(inm)*kprimes(n+1) + bs(inm)*iprimes(n+1);
			
			ene = ene + c1*c2* 4*pi/(2*n+1);
		end
	end
	ene = ene * ra^2 * 0.5;
	
	eint = (ene+(kappa+1)*2*pi)*2;
	
	disp(['h=',num2str(norm(pa-pb)-ra-rb)]);
	disp(['kh=',num2str(kappa*(norm(pa-pb)-ra-rb))]);
	disp(['ene=',num2str(ene)])
	disp(['eint=',num2str(eint)])
end

if 1
	% interaction force
	
	% for this sphere
	cs = zeros(size(coefa));
	for n = 0:nmax
	for m = -n:n
		inm = sh_sub2ind(n,m);
		cs(inm) = as(inm)*kprimes(n+1) + bs(inm)*iprimes(n+1);
	end
	end
	
	% can be shared by all
	alphas = zeros(npole);
	betas = zeros(npole);
	Ns = zeros(npole);
	for n = 0:nmax
	for m = -n:n
		inm = sh_sub2ind(n,m);
		
		absm = abs(m);
		
		%
		if m >= 0
			alphas(inm) = 1.0;
		else
			alphas(inm) = (-1)^absm * factorial(n-absm) / factorial(n+absm);
		end
		
		%
		betas(inm) = 2 * factorial(n+absm) / (2*n+1) / factorial(n-absm);
		
		% 
		Ns(inm) = (-1)^m * sqrt(factorial(n-absm) / factorial(n+absm));
		
	end
	end
	
	fint = [0,0,0]';
	for n = 0:nmax
	for m = -n:n
		inm = sh_sub2ind(n,m);
		inm1 = sh_sub2ind(n,-m);
		
		betanm = betas(inm);
		
		for p = 0:nmax
		for q = -p:p
			ipq = sh_sub2ind(p,q);
			
			% x,y force
			if q==-m-1 || q==-m+1
				
				c_sin_theta = 0;
				if     p==n-1 && q==-m-1
					c_sin_theta = -alphas(inm1)/alphas(ipq) * 1/(2*n-1) * betanm;
				elseif p==n+1 && q==-m-1
					c_sin_theta = alphas(inm1)/alphas(ipq) * 1/(2*n+3) * betanm;
				elseif p==n-1 && q==-m+1
					c_sin_theta = alphas(inm1)/alphas(ipq) * (n+m-1)*(n+m)/(2*n-1) * betanm;
				elseif p==n+1 && q==-m+1
					c_sin_theta = -alphas(inm1)/alphas(ipq) * (n-m+1)*(n-m+2)/(2*n+3) * betanm;
				end
				
				c_cos_phi = 0;
				c_sin_phi = 0;
				if     q == -m-1
					c_cos_phi = pi;
					c_sin_phi = -pi * 1i;
				elseif q == -m+1
					c_cos_phi = pi;
					c_sin_phi = pi * 1i;
				end
				
				cnm = cs(inm);
				cpq = cs(ipq);
				Nnm = Ns(inm);
				Npq = Ns(ipq);
				
				fint(1) = fint(1) + cnm*cpq*Nnm*Npq * c_sin_theta * c_cos_phi;
				fint(2) = fint(2) + cnm*cpq*Nnm*Npq * c_sin_theta * c_sin_phi;
			end
			
			% z force
			if q==-m
				
				c_cos_theta = 0;
				if     p==n-1
					c_cos_theta = alphas(inm1)/alphas(ipq) * (n+m)/(2*n-1) * betanm;
				elseif p==n+1
					c_cos_theta = alphas(inm1)/alphas(ipq) * (n-m+1)/(2*n+3) * betanm;
				end
				
				c_phi = pi*2;
				
				fint(3) = fint(3) + cs(inm)*cs(ipq)*Ns(inm)*Ns(ipq) * c_cos_theta * c_phi;
			end
		end
		end
	end
	end
	fint = fint .* (0.5*ra^2);
	
	(fint)
	fnorm = norm(fint)
	fint ./ fnorm
	(pa-pb) ./ norm(pa-pb)
end
























