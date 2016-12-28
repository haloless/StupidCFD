
function [s1,s1a,s1b] = EwaldS1(x,y,z,lx,ly,lz,alpha,nmax,kmax)

twoPiI = pi * 2.0i;
tau0 = lx * ly * lz;
coef = pi*alpha / gamma(1);

s1a = 0;
s1b = 0;

for ii = -nmax:nmax
for jj = -nmax:nmax
for kk = -nmax:nmax
	dx = x - ii*lx;
	dy = y - jj*ly;
	dz = z - kk*lz;
	dr2 = dx^2 + dy^2 + dz^2;
	
	ss = FuncPhiNhalf(pi*dr2/alpha);
	s1a = s1a + ss;
end
end
end

s1a = coef * (s1a * tau0 * alpha^(-1.5) - 1);

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
	ff = FuncPhi0(pi*alpha*k2);
	s1b = s1b + ee*ff;
end
end
end

s1b = real(s1b) * coef;

s1 = s1a+s1b;

s1 = s1 / (pi*tau0);

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

% function [val] = FuncPhiNhalf(x)

% return
% end

