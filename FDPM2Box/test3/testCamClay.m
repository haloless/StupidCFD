clear;

bm1 = [1;1;1;0;0;0];

par = struct();

par.E = 100;
par.nu = 0.1;

par.ptens = 0.1;
par.beta = 0.25;
par.M = 5;
par.Hslope = 1.2;
% par.Hslope = 0;

%
par.a0 = par.ptens / (1+par.beta);
% par.a0 = 3;

%
epsEP = zeros(6,1);
epsE = zeros(6,1);
sigma = zeros(6,1);
alpha = 0;


if 1
	hfig = figure;
	
	pp = -0.1:0.05:0.1; qq = -par.M*(pp-par.ptens);
	plot(pp,qq,'-b', pp,zeros(size(pp)),'-k', zeros(size(qq)),qq,'-k');
	xlabel('p'); ylabel('q');
	
	hold on;
	aa = par.a0 + alpha*par.Hslope;
	pa = par.ptens - aa;
	tt = 0:0.01:pi/2; xx1 = aa*cos(tt); yy1 = par.M*aa*sin(tt);
	tt = tt + pi/2; xx2 = par.beta*aa*cos(tt); yy2 = par.M*aa*sin(tt);
	plot(pa+xx1,yy1,'g-',pa+xx2,yy2,'r-');
	hold off;
	% return;
end

%
epsconsmax = -0.02;

nstep = 10;
for istep = 1:nstep
	epsEold = epsE;
	alphaold = alpha;
	
	deps = zeros(6,1);
	deps(1:3) = epsconsmax / nstep;
	
	epsEtr = epsEold + deps;
	[dgam,sigma,epsE,alpha,Dalg] = material3dCamClay(par, epsEtr,alphaold);
	
	if 1
		figure(hfig);
		hold on;
		
		p = mean(sigma(1:3));
		s = sigma - p*bm1;
		q = sqrt(3/2) * sqrt(sum(s(1:3).^2) + 2*sum(s(4:6).^2));
		plot(p,q,'x');
		
		aa = par.a0 + alpha*par.Hslope;
		pa = par.ptens - aa;
		tt = 0:0.01:pi/2; xx1 = aa*cos(tt); yy1 = par.M*aa*sin(tt);
		tt = tt + pi/2; xx2 = par.beta*aa*cos(tt); yy2 = par.M*aa*sin(tt);
		plot(pa+xx1,yy1,'g-',pa+xx2,yy2,'r-');
		
		hold off;
		pause;
	end
	if 1
		Ddiff = zeros(6);
		dh = 1.0e-6;
		for dir = 1:6
			epsEtrp = epsEtr; epsEtrp(dir) = epsEtrp(dir) + dh;
			[~,sigmap] = material3dCamClay(par,epsEtrp,alphaold);
			epsEtrm = epsEtr; epsEtrm(dir) = epsEtrm(dir) - dh;
			[~,sigmam] = material3dCamClay(par,epsEtrm,alphaold);
			Ddiff(:,dir) = (sigmap - sigmam) ./ (dh*2);
		end
		Dalg
		Ddiff
		pause;
	end
end


%
epsxxmax = -0.1;
nstep = 50;
depsxx = (epsxxmax - epsconsmax) / nstep;

for istep = 1:nstep
	epsEold = epsE;
	alphaold = alpha;
	
	deps = zeros(6,1);
	deps(1) = depsxx;
	deps(2:3) = -depsxx*0.8;
	
	epsEtr = epsEold + deps;
	[dgam,sigma,epsE,alpha,Dalg] = material3dCamClay(par, epsEtr,alphaold);
	
	if 1
		figure(hfig);
		hold on;
		
		p = mean(sigma(1:3));
		s = sigma - p*bm1;
		q = sqrt(3/2) * sqrt(sum(s(1:3).^2) + 2*sum(s(4:6).^2));
		plot(p,q,'x');
		
		aa = par.a0 + alpha*par.Hslope;
		pa = par.ptens - aa;
		tt = 0:0.01:pi/2; xx1 = aa*cos(tt); yy1 = par.M*aa*sin(tt);
		tt = tt + pi/2; xx2 = par.beta*aa*cos(tt); yy2 = par.M*aa*sin(tt);
		plot(pa+xx1,yy1,'g-',pa+xx2,yy2,'r-');
		
		hold off;
		pause;
	end
	if 1
		Ddiff = zeros(6);
		dh = 1.0e-8;
		for dir = 1:6
            % NOTE if approaching the critical point, should use one-sided finite difference
            if 0
                epsEtrp = epsEtr; epsEtrp(dir) = epsEtrp(dir) + dh;
                [~,sigmap] = material3dCamClay(par,epsEtrp,alphaold);
                epsEtrm = epsEtr; epsEtrm(dir) = epsEtrm(dir) - dh;
                [~,sigmam] = material3dCamClay(par,epsEtrm,alphaold);
                Ddiff(:,dir) = (sigmap - sigmam) ./ (dh*2);
            else
                aa = par.a0 + alpha*par.Hslope;
                pa = par.ptens - aa;
                if p>=pa
                    epsEtrp = epsEtr; epsEtrp(dir) = epsEtrp(dir) + dh;
                    [~,sigmap] = material3dCamClay(par,epsEtrp,alphaold);
                    Ddiff(:,dir) = (sigmap-sigma) ./ dh;
                else
                    epsEtrm = epsEtr; epsEtrm(dir) = epsEtrm(dir) - dh;
                    [~,sigmam] = material3dCamClay(par,epsEtrm,alphaold);
                    Ddiff(:,dir) = (sigma-sigmam) ./ dh;
                end
            end
		end
		Dalg
		Ddiff
		pause;
	end
end
