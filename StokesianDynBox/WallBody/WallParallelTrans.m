%% Resistance for sphere-wall translation in parallel direction
%% Solve the tri-diagonal system of (O'Neill, 1964)
%%

function [fstar,gstar] = WallParallelTrans(ha,nmax)

% truncation of series
if ~exist('nmax','var')
	nmax = 1000;
end
nmax2 = nmax + 2;

alpha = acosh(ha);
sinha = sinh(alpha);
cotha = coth(alpha);
cosecha = 1 / sinha;

% generate k_n
ks = zeros(nmax2,1);
for n = 0:nmax2-1
    ks(n+1) = kfunc(n,alpha);
end

% solve A_n
mat = sparse(nmax,nmax);
rhs = zeros(nmax,1);
for i = 1:nmax
    ca = (2*i-1)*ks(i) - (2*i-3)*ks(i+1);
    cb = (2*i+5)*ks(i+1) - (2*i+3)*ks(i+2);
    
    c0 = ca * (i-1)/(2*i-1);
    c1 = -ca*i/(2*i+1) - cb*(i+1)/(2*i+1);
    c2 = cb * (i+2)/(2*i+3);
    
    if i-1 >= 1
        mat(i,i-1) = c0;
    end
    mat(i,i) = c1;
    if i+1 <= nmax  
        mat(i,i+1) = c2;
    end
    
    a1 = (i+0.5) * alpha;
    a2 = (i-0.5) * alpha;
    a3 = (i+1.5) * alpha;
    rhs(i) = sqrt(2) * (2*coth(a1) - coth(a2) - coth(a3));
end

As = mat \ rhs;

% F*
fstar = 0;
for n = 0:nmax-1
    En = Efunc(n,alpha,ks,As);
    Cn = Cfunc(n,alpha,ks,As);
    
    fstar = fstar + En + n*(n+1)*Cn;
end
fstar = 1/6*sqrt(2) * sinh(alpha) * fstar;

% G*
gstar = 0;
for n = 0:nmax-1
	
    Bn = Bfunc(n,alpha,ks,As);
    Cn = Cfunc(n,alpha,ks,As);
    Dn = Dfunc(n,alpha,ks,As);
    En = Efunc(n,alpha,ks,As);
	
	An = 0;
	if n > 0
		An = As(n);
	end
	
	ea = exp(-(2*n+1)*alpha);
	gstar = gstar + (2+ea) * (n*(n+1)*(2*An+Cn*cotha) - (2*n+1-cotha)*En);
	gstar = gstar + (2-ea) * (n*(n+1)*Bn*cotha - (2*n+1-cotha)*Dn);
end
gstar = gstar / (12*sqrt(2)*cosecha^2);
% we normalize by 6*pi*mu*a^3
gstar = gstar * 4/3;



return
end


function [kn] = kfunc(n,alpha)
	kn = (n+0.5) * coth((n+0.5)*alpha) - coth(alpha);
return
end

function [En] = Efunc(n,alpha,ks,As)
	En = 2*sqrt(2)*exp(-(n+0.5)*alpha) / sinh((n+0.5)*alpha);
    if n >= 2
        En = En + ks(n+1) * (n-1)*n/(2*n-1) * As(n-1);
    end
    En = En - ks(n+1) * (n+1)*(n+2)/(2*n+3) * As(n+1);
	% En = 0;
	% if n >= 2
		% En = 0.5 * (As(n-1)-As(n+1));
	% end
return
end

function [Cn] = Cfunc(n,alpha,ks,As)
	Cn = 0;
    if n >= 1
        if n >= 2
            Cn = Cn + (n-1)/(2*n-1)*As(n-1);
        end
        Cn = Cn - As(n) + (n+2)/(2*n+3)*As(n+1);
        Cn = -2*ks(n+1) * Cn;
    end
return
end

function [Bn] = Bfunc(n,alpha,ks,As)
	Bn = 0;
	if n >= 1
		if n >= 2
			Bn = Bn + (n-1)*As(n-1);
		end
		Bn = Bn - (2*n+1)*As(n) + (n+2)*As(n+1);
	end
return
end

function [Dn] = Dfunc(n,alpha,ks,As)
	Dn = 0;
	if n >= 2
		Dn = Dn - 0.5*(n-1)*n*As(n-1);
	end
	Dn = Dn + 0.5*(n+1)*(n+2)*As(n+1);
return
end













