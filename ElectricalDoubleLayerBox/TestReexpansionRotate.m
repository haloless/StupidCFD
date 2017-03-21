
clear;

ImUnit = 1i;

kappa = 1.0;
% kappa = 4.0;

func_reg = @(n,m,r,theta,phi) r^n * AdaptModSphBesselI(n,kappa*r) * SphHarmY(n,m,theta,phi);
func_sing = @(n,m,r,theta,phi) exp(-kappa*r) / r^(n+1) * AdaptModSphBesselK(n,kappa*r) * SphHarmY(n,m,theta,phi);


pa = [ 0; 0; 0 ];
% pb = [ 0; 0; 4 ];
% pb = [ 1; 0; 4 ];
% pb = [ 0; 1; 4 ];
pb = [ 1; 1; 4 ];
pb = -pb;
% pb = rand(3,1)

% ptmp = pa;
% pa = pb;
% pb = ptmp;

pab = pa - pb;
% pab = pb - pa;

% create original->coaxial rotation
ReExpSimpleRot;
% return

%
[r,theta,phi] = sh_cart2sph(pab(1),pab(2),pab(3));
kr = kappa * r;

if 1
    theta = thetaprime;
    phi = phiprime;
    % phi = phiprime - pi/2;
    % phi = phi - pi/2;
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
ptest = ptest ./ norm(ptest);
ptest = ptest + pa;

disp(['ptest=',num2str(ptest')]);

[ra,thetaa,phia] = sh_cart2sph(ptest(1)-pa(1),ptest(2)-pa(2),ptest(3)-pa(3));
[rb,thetab,phib] = sh_cart2sph(ptest(1)-pb(1),ptest(2)-pb(2),ptest(3)-pb(3));


% figure;
% plot3([pa(1),pb(1)],[pa(2),pb(2)],[pa(3),pb(3)],'x-')

% nmax = 2;
nmax = 4;
% nmax = 8;
% nmax = 10;
% nmax = 12;
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


if 0
    disp(['current -> coaxial']);
    phat = Qmat * (ptest-pa);
    [rhat,thetahat,phihat] = sh_cart2sph(phat(1),phat(2),phat(3));
    
    ytest = zeros(npole,1);
    yhat = zeros(npole,1);
    
    for n = 0:nmax
        for m = -n:n
            ytest(sh_sub2ind(n,m)) = SphHarmY(n,m,thetaa,phia);
        end
    end
    
    for n = 0:nmax
        for m = -n:n
            ynm = 0;
            for s = -n:n
                Rnms = Rmat(sh_sub2ind(n,m),sh_sub2ind(n,s));
                yns = Rnms * SphHarmY(n,s,thetahat,phihat);
                ynm = ynm + yns;
            end
            yhat(sh_sub2ind(n,m)) = ynm;
        end
    end
    
    
    return
end

if 0
    disp(['coaxial -> current']);
    phat = Qmat * (ptest-pa);
    [rhat,thetahat,phihat] = sh_cart2sph(phat(1),phat(2),phat(3));
    
    ytest = zeros(npole,1);
    yhat = zeros(npole,1);
    
    for n = 0:nmax
        for m = -n:n
            yhat(sh_sub2ind(n,m)) = SphHarmY(n,m,thetahat,phihat);
        end
    end
    
    Rhatmat = Rmat';
    
    for n = 0:nmax
        for m = -n:n
            ynm = 0;
            for s = -n:n
                Rnms = Rhatmat(sh_sub2ind(n,m),sh_sub2ind(n,s));
                yns = Rnms * SphHarmY(n,s,thetaa,phia);
                ynm = ynm + yns;
            end
            ytest(sh_sub2ind(n,m)) = ynm;
        end
    end
    
    
    return
end

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

% bcoef = rand(npole,1);
% bcoef = (1.0 ./ (1:npole).^2)';
% bcoef = exp(-(1:npole).^2)';
bcoef = zeros(npole,1);
for n = 0:nmax
for m = -n:n
	inm = sh_sub2ind(n,m);
	% bcoef(inm) = exp(-n) * (m^2+1);
	bcoef(inm) = 1 / bkhat(n+1);
    % bcoef(inm) = rand(1);
end
end

%
Tmat = Rmat' * Smat * Rmat;
% Tmat = Rmat * Smat * Rmat';
% Tmat = Tmat';
% Tmat = Smat;
lcoef = Tmat * bcoef;


% bval = 0;
% for n = 0:nmax
	% for m = -n:n
		% inm = sh_sub2ind(n,m);
		% valnm = bkhat(n+1) * bcoef(inm) * bshy(inm);
		% bval = bval + valnm;
	% end
% end

% aval = 0;
% for n = 0:nmax
	% for m = -n:n
		% inm = sh_sub2ind(n,m);
		% valnm = aihat(n+1) * lcoef(inm) * ashy(inm);
		% aval = aval + valnm;
	% end
% end

% disp(['bval=',num2str(bval)]);
% disp(['aval=',num2str(aval)]);

bval = zeros(npole,1);
for n = 0:nmax
for m = -n:n
	inm = sh_sub2ind(n,m);
	bval(inm) = bcoef(inm) * bkhat(n+1) * bshy(inm);
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
		
		valnm = bcoef(ipq) * Tmat(inm,ipq) * aihat(n+1) * ashy(inm);
		val = val + valnm;
	end
	end
	
	aval(ipq) = val;
end
end

[bval(1:4), aval(1:4)]

