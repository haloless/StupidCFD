function [f,df] = ReExpFuncRegular(nmax,r,kappa, add_adapt_coef)
% Regular solution
% ihat(r;kappa)
% dihat/dr
% 

if ~exist('add_adapt_coef','var')
	add_adapt_coef = 1;
end


% number of terms
nmax1 = nmax + 1;

% actual variable 
x = r * kappa;
% cutoff value
xsmall = 1.0e-8;
% x = max(x,xsmall);

% pre-calculate values
sinhx = sinh(x);
coshx = cosh(x);

% 
f = zeros(nmax1,1);
df = zeros(nmax1,1);
for n = 0:nmax
	if n == 0
		% initial value for n=0,1
		fm = 0;
		if x > xsmall
			fn = sinhx / x;
			fp = (x*coshx-sinhx) / x^2;
		else
			% for x~0
			fn = 1;
			fp = 0;
		end
		
		% derivative i0' equal to i1
		dfn = fp;
	else
		if x > xsmall
			fp = fm - (2*n+1)/x*fn;
		else
			fp = 0;
		end
		
		dfn = n/(2*n+1)*fm + (n+1)/(2*n+1)*fp;
	end
	
	f(n+1) = fn;
	df(n+1) = dfn;
	
	fm = fn;
	fn = fp;
end

if add_adapt_coef
	%
	for n = 0:nmax
		coef = double_factorial(2*n+1) / kappa^n;
		f(n+1) = f(n+1) * coef;
		df(n+1) = df(n+1) * coef*kappa;
	end
end


return
end



