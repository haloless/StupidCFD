
clear;

ImUnit = 1i;

% kappa = 1.0;
kappa = 4.0;

func_reg = @(n,m,r,theta,phi) r^n * AdaptModSphBesselI(n,kappa*r) * SphHarmY(n,m,theta,phi);
func_sing = @(n,m,r,theta,phi) exp(-kappa*r) / r^(n+1) * AdaptModSphBesselK(n,kappa*r) * SphHarmY(n,m,theta,phi);


pa = [ 0; 0; 0 ];
% pb = [ 0; 0; 4 ];
% pb = [ 1; 0; 4 ];
% pb = [ 0; 1; 4 ];
pb = [ -1; -1; 4 ];


ptest = rand(3,1);
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
% reexpansion matrix
%
[Tmat,Rmat,Smat] = ReExpMat(pb,pa, kappa,nmax);


%
% test
%
bcoef = zeros(npole,1);
for n = 0:nmax
for m = -n:n
	inm = sh_sub2ind(n,m);
	bkhat = exp(-kappa*rb)/rb^(n+1) * AdaptModSphBesselK(n,kappa*rb);
	% bcoef(inm) = rand(1) / bkhat;
	bcoef(inm) = (rand(1)+ImUnit*rand(1))*100 / bkhat;
end
end


%
areg = zeros(npole,1);
breg = zeros(npole,1);
asing = zeros(npole,1);
bsing = zeros(npole,1);
for n = 0:nmax
for m = -n:n
    inm = sh_sub2ind(n,m);
    areg(inm) = func_reg(n,m,ra,thetaa,phia);
	breg(inm) = func_reg(n,m,rb,thetab,phib);
	asing(inm) = func_sing(n,m,ra,thetaa,phia);
    bsing(inm) = func_sing(n,m,rb,thetab,phib);
end
end

bval = bcoef' * bsing;
disp(['bval=',num2str(bval)]);

acoef = Tmat * bcoef;
aval = acoef' * areg;
disp(['aval=',num2str(aval)]);


%
if 1
	% compare expanded base functions
	areexp = Tmat' * areg;
	[bsing(1:4), areexp(1:4)]
	
	% breexp = -Tmat * breg;
	% [asing(1:4), breexp(1:4)]
end

if 0
	% compare reverse expansion
	Trev = ReExpMat(pa,pb,kappa,nmax);
	breexp = Trev' * breg;
	[asing(1:4), breexp(1:4)]
end






