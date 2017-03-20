
clear all;

ImUnit = 1i;

kappa = 1.0;
% kappa = 4.0;

func_reg = @(n,m,r,theta,phi) r^n * AdaptModSphBesselI(n,kappa*r) * SphHarmY(n,m,theta,phi);
func_sing = @(n,m,r,theta,phi) exp(-kappa*r) / r^(n+1) * AdaptModSphBesselK(n,kappa*r) * SphHarmY(n,m,theta,phi);


pa = [ 0; 0; 0 ];
% pb = [ 0; 0; 4 ];
% pb = [ 1; 0; 4 ];
% pb = [ 0; 1; 4 ];
pb = [ -1; -1; 4 ];
% pb = rand(3,1);
% pb = pb ./ norm(pb) * 5;
% pb = -pb;

pab = pa - pb;

% create original->coaxial rotation
ReExpSimpleRot;
% return

%
[r,theta,phi] = sh_cart2sph(pab(1),pab(2),pab(3));
kr = kappa * r;

if 1
    theta = thetaprime;
    phi = phiprime;
end






% ptest = [ 0.1; 0.1; 0.1 ];
% ptest = [ 2; 0.1; 0.5 ];
% ptest = [ 0.2; 0.1; 0.5 ];
% ptest = [ 1; 0.5; 0.2 ];
% ptest = [ 1; 1; 2.5 ];
% ptest = [ 1; 0; 2.5 ];
% ptest = [ 2; 2; 1.5 ];
% ptest = [ 2; 3; 1.5 ];
% ptest = [ 4; 3; 2 ];
ptest = rand(3,1);
% ptest = -ptest;
ptest = ptest ./ norm(ptest);
ptest = ptest + pa;

disp(['ptest=',num2str(ptest')]);

[ra,thetaa,phia] = sh_cart2sph(ptest(1)-pa(1),ptest(2)-pa(2),ptest(3)-pa(3));
[rb,thetab,phib] = sh_cart2sph(ptest(1)-pb(1),ptest(2)-pb(2),ptest(3)-pb(3));

% nmax = 2;
% nmax = 4;
% nmax = 8;
% nmax = 10;
nmax = 12;
% nmax = 16;
% nmax = 24;
% nmax = 32;
nmax1 = nmax + 1;
npole = nmax1^2;
disp(['nmax=',int2str(nmax), '; npole=',int2str(npole)]);


%
% rotation matrix
%
ReExpCalcR;


%
% scale matrix
%
ReExpCalcS;


ashy = zeros(npole,1);
bshy = zeros(npole,1);
for n = 0:nmax
for m = -n:n
	inm = sh_sub2ind(n,m);
	ashy(inm) = SphHarmY(n,m,thetaa,phia);
	bshy(inm) = SphHarmY(n,m,thetab,phib);
end
end

aihat = zeros(nmax1,1);
bkhat = zeros(nmax1,1);
for n = 0:nmax
	bkhat(n+1) = exp(-kappa*rb)/rb^(n+1) * AdaptModSphBesselK(n,kappa*rb);
	aihat(n+1) = ra^n * AdaptModSphBesselI(n,kappa*ra);
end

bcoef = zeros(npole,1);
for n = 0:nmax
for m = -n:n
	inm = sh_sub2ind(n,m);
	% bcoef(inm) = 1 / bkhat(n+1);
	bcoef(inm) = rand(1) / bkhat(n+1);
end
end


%
areg = zeros(npole,1);
bsing = zeros(npole,1);
for n = 0:nmax
for m = -n:n
    inm = sh_sub2ind(n,m);
    areg(inm) = func_reg(n,m,ra,thetaa,phia);
    bsing(inm) = func_sing(n,m,rb,thetab,phib);
end
end

bval = bcoef' * bsing;
disp(['bval=',num2str(bval)]);

% 1st rotate
[rahat,thetaahat,phiahat] = sh_cart2sph(Qmat*(ptest-pa));
[rbhat,thetabhat,phibhat] = sh_cart2sph(Qmat*(ptest-pb));

bsinghat = zeros(npole,1);
for n = 0:nmax
for m = -n:n
    bsinghat(sh_sub2ind(n,m)) = func_sing(n,m,rbhat,thetabhat,phibhat);
end
end
bcoefhat = Rmat' * bcoef;
bvalhat = bcoefhat' * bsinghat;
disp(['bvalhat=',num2str(bvalhat)]);

% 2nd translate
areghat = zeros(npole,1);
for n = 0:nmax
for m = -n:n
    areghat(sh_sub2ind(n,m)) = func_reg(n,m,rahat,thetaahat,phiahat);
end
end
acoefhat = Smat * bcoefhat;
avalhat = acoefhat' * areghat;
disp(['avalhat=',num2str(avalhat)]);

% 3rd rotate
acoef = Rmat * acoefhat;
aval = acoef' * areg;
disp(['aval=',num2str(aval)]);

%
Tmat = Rmat * Smat * Rmat';
acoef = Tmat * bcoef;
aval = acoef' * areg;
disp(['aval=',num2str(aval)]);

% return

%
areexp = Tmat' * areg;
% return

bval = zeros(npole,1);
for n = 0:nmax
for m = -n:n
	inm = sh_sub2ind(n,m);
	% bval(inm) = bcoef(inm) * bkhat(n+1) * bshy(inm);
	% bval(inm) = bcoef(inm) * bsing(inm);
	bval(inm) = bsing(inm);
end
end

aval = zeros(npole,1);
for p = 0:nmax
for q = -p:p
	ipq = sh_sub2ind(p,q);
	
	val = 0;
	for n = 0:nmax
	for m = -n:n
		inm = sh_sub2ind(n,m);
		
		% valnm = bcoef(ipq) * Tmat(inm,ipq) * aihat(n+1) * ashy(inm);
		% valnm = bcoef(ipq) * Tmat(inm,ipq) * areg(inm);
		valnm = conj(Tmat(inm,ipq)) * areg(inm);
		% valnm = Tmat(inm,ipq) * areg(inm);
		val = val + valnm;
	end
	end
	
	aval(ipq) = val;
end
end

[bval(1:4), aval(1:4)]

