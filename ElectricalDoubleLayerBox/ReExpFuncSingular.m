function [f,df] = ReExpFuncSingular(nmax,r,kappa)
% Singular solution 
% khat(r;kappa) = 
% dkhat/dr = 
% This function is OK for kappa=0
%

% number of terms
nmax1 = nmax + 1;

% actual variable 
x = r * kappa;

% 
f = zeros(nmax1,1);
df = zeros(nmax1,1);
for n = 0:nmax
	if n == 0
		% initial value for n=0,1
		fm = 0;
		fn = exp(-x) ./ r;
		fp = exp(-x) ./ r^2 .* (x+1);
		
		% derivative n=0 is -k1
		dfn = -fp;
	else
		% next function
		fp = kappa^2 * fm + (2*n+1)/r * fn;
		
		% this derivative
		dfn = -n/(2*n+1)*kappa^2*fm - (n+1)/(2*n+1)*fp;
	end
	
	f(n+1) = fn;
	df(n+1) = dfn;
	
	fm = fn;
	fn = fp;
end

if 1
	% add (2n-1)!! denominator
	for n = 0:nmax
		coef = double_factorial(2*n-1);
		f(n+1,:) = f(n+1,:) ./ coef;
		df(n+1,:) = df(n+1,:) ./ coef;
	end
end



return
end


