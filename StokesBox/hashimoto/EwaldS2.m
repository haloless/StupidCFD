
function [s2,s2a,s2b] = EwaldS2(x,y,z,lx,ly,lz,alpha,nmax,kmax)

twoPiI = pi * 2.0i;
tau0 = lx * ly * lz;
coef = pi^2 * alpha^2 / gamma(2);

s2a = 0;
s2b = 0;

for ii = -nmax:nmax
for jj = -nmax:nmax
for kk = -nmax:nmax
	dx = x - ii*lx;
	dy = y - jj*ly;
	dz = z - kk*lz;
	dr2 = dx^2 + dy^2 + dz^2;
	
	xx = pi*dr2/alpha;
	ff = 2*exp(-xx) - 2*xx*FuncPhiNhalf(xx);
	s2a = s2a + ff;
end
end
end

s2a = coef * (s2a * tau0 * alpha^(-1.5) - 0.5);

for ii = -kmax:kmax
for jj = -kmax:kmax
for kk = -kmax:kmax
	if (ii==0 && jj==0 && kk==0)
		continue;
	end
	
	kx = ii / lx;
	ky = jj / ly;
	kz = kk / lz;
	k2 = kx^2 + ky^2 + kz^2;
	kr = kx*x + ky*y + kz*z;
	
	ee = exp(twoPiI*kr);
	xx = pi*alpha*k2;
	ff = (exp(-xx)+FuncPhi0(xx))/xx;
	s2b = s2b + ee*ff;
end
end
end

s2b = real(s2b) * coef;

s2 = s2a+s2b;

s2 = -s2 / (4*pi^3*tau0);

return
end

function [val] = FuncPhi0(x)
val = exp(-x) ./ x;
return
end

function [val] = FuncPhiNhalf(x)
sqrtx = sqrt(x);
val = sqrt(pi) ./ sqrtx .* erfc(sqrtx);
return
end

