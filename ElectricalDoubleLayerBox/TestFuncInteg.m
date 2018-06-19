
function [] = TestFuncInteg()

if 0
	%
	% test integrate theta-part with sin(theta)
	%
	
	fun = @(n,m,p,q,x) AssocLegendrePhat(n,m,x) .* AssocLegendrePhat(p,q,x) .* sqrt(1-x.^2)
	
	nmax = 5;
	
	for n = 0:nmax
	for m = -n:n
		
		betanm = PIntCoef(n,abs(m));
		
		for p = 0:nmax
		for q = -p:p
			if q==-m-1 || q==-m+1
				uint = integral(@(x) fun(n,m,p,q,x), -1.0,1.0);
				
				uana = 0;
				if     p==n-1 && q==-m-1
					uana = -PGenCoef(n,-m)/PGenCoef(p,q) * 1/(2*n-1) * betanm;
				elseif p==n+1 && q==-m-1
					uana = PGenCoef(n,-m)/PGenCoef(p,q) * 1/(2*n+3) * betanm;
				elseif p==n-1 && q==-m+1
					uana = PGenCoef(n,-m)/PGenCoef(p,q) * (n+m-1)*(n+m)/(2*n-1) * betanm;
				elseif p==n+1 && q==-m+1
					uana = -PGenCoef(n,-m)/PGenCoef(p,q) * (n-m+1)*(n-m+2)/(2*n+3) * betanm;
				end
				
				uerr = abs(uana-uint);
				ok = uerr < 1.0e-10;
				assert(ok);
				
				fprintf('n=%d,m=%d,p=%d,q=%d,ok=%d,uint=%e,uana=%e\n', n,m,p,q, ok, uint,uana);
			end
		end
		end
	end
	end
		
end

if 1
	%
	% test integrate theta-part with cos(theta)
	%
	
	fun = @(n,m,p,q,x) AssocLegendrePhat(n,m,x) .* AssocLegendrePhat(p,q,x) .* x
	
	nmax = 5;
	
	for n = 0:nmax
	for m = -n:n
		
		betanm = PIntCoef(n,abs(m));
		
		for p = 0:nmax
		for q = -p:p
			if q==-m
				uint = integral(@(x) fun(n,m,p,q,x), -1.0,1.0);
				
				uana = 0;
				if     p==n-1 && q==-m
					uana = PGenCoef(n,-m)/PGenCoef(p,q) * (n+m)/(2*n-1) * betanm;
				elseif p==n+1 && q==-m
					uana = PGenCoef(n,-m)/PGenCoef(p,q) * (n-m+1)/(2*n+3) * betanm;
				end
				
				uerr = abs(uana-uint);
				ok = uerr < 1.0e-10;
				assert(ok);
				
				fprintf('n=%d,m=%d,p=%d,q=%d,ok=%d,uint=%e,uana=%e\n', n,m,p,q, ok, uint,uana);
			end
		end
		end
	end
	end
		
end

return
end

% function [] = Tes

function [pnm] = AssocLegendrePhat(n,m,x)
	absm = abs(m);
	pnm = AssocLegendreP(n,absm,x);
	return
end

function [pnm] = AssocLegendreP(n,m,x)
	% assert(-n<=m && m<=n)
	% assert(-1<=x && x<=1)
	
	am = abs(m);
	
	szx = size(x);
	
	pnm = legendre(n,x);
	pnm = pnm(am+1,:)';
	
	if 1
		pnm = pnm * PGenCoef(n,m);
	end
	
	pnm = reshape(pnm,szx);
	
	return
end

function [coef] = PIntCoef(n,m)
	coef = 2 * factorial(n+m) / (2*n+1) / factorial(n-m);
	return
end

function [anm] = PGenCoef(n,m)
	if m >= 0
		anm = 1.0;
	else
		m = -m;
		anm = (-1)^m * factorial(n-m) / factorial(n+m);
	end
	return
end
