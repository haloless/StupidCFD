
clear all;

a = 1.0;
% 
kappa = 2.0 / a;
% kappa = 1.0 / a;
% kappa = 0.5 / a;
% kappa = 0.1 / a;
% kappa = 0.01 / a;
% 
ak = a * kappa;
lambda = 1 / kappa;


% separation
% H = 1.0 * a;
% H = 0.2 * a;
H = 0.1 * a;
% H = 0.01 * a;
% center distance
R = H + a + a;


cutoff = 20;
cutoff2 = 15;


Kak = zeros(cutoff,1);
Iak = zeros(cutoff,1);
for n = 1:cutoff
	nn = n - 1;
	Kak(n) = ModSphBesselK(nn, ak);
	Iak(n) = ModSphBesselI(nn, ak);
end

Bmat = zeros(cutoff,cutoff);
for n = 1:cutoff
for m = 1:cutoff
	nn = n - 1;
	mm = m - 1;
	
	Bnm = 0;
	numax = min(nn,mm);
	for nu = 0:numax
		Bnm = Bnm + FuncA(nu,nn,mm) * ModSphBesselK(nn+mm-2*nu,ak*R);
	end
	
	Bmat(n,m) = Bnm;
end
end

Lmat = zeros(cutoff,cutoff);
for j = 1:cutoff
for n = 1:cutoff
	jj = j - 1;
	nn = n - 1;
	
	Ljn = (2*jj+1) * Bmat(n,j) * Iak(j) / Kak(n);
	
	Lmat(j,n) = Ljn;
end
end

Imat = eye(cutoff);
mat = Lmat + Imat;

rhs = zeros(cutoff,1);
rhs(1) = 1;

sol = mat \ rhs;

acoef = sol ./ Kak;


if 0
    [xs,ys] = ndgrid(-a*2:0.1:R+a*2,0:0.1:a*2);
    nx = size(xs,1);
    ny = size(xs,2);
    rs = sqrt(xs.^2 + ys.^2);
    ts = atan2(ys,xs);

    ok = ones(nx,ny);
    ok((xs-0).^2+(ys-0).^2<a^2) = 0;
    ok((xs-R).^2+(ys-0).^2<a^2) = 0;

    phis = zeros(nx,ny);

    for i = 1:nx
    for j = 1:ny
        if ok == 0 
            continue;
        end
        
        xx = xs(i,j);
        yy = ys(i,j);
        rr = rs(i,j);
        rk = rr * kappa;
        theta = ts(i,j);
        ct = cos(theta);
        
        uu = 0;
        for n = 1:cutoff2
            nn = n - 1;
            
            an = acoef(n);
            un = 0;
            for m = 1:cutoff2
                mm = m - 1;
                un = un + an * (2*mm+1) * Bmat(n,m) * ModSphBesselI(mm,rk) * LegendrePoly(mm,ct);
            end
            
            un = un + an * ModSphBesselK(nn,rk) * LegendrePoly(nn,ct);
            
            uu = uu + un;
        end
        
        phis(i,j) = uu;
    end
    end

    if 1
    phis(ok==0) = nan;
    figure; contourf(xs,ys,phis,0:0.05:1); axis equal
    end

end