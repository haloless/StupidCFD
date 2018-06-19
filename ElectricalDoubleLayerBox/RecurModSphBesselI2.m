function [I] = RecurModBesselI2(nmax,x)

I = zeros(nmax+1,1);

sinhx = sinh(x);
coshx = cosh(x);

for n = 0:nmax
	if n == 0
		% initial value for n=0,1
		im = 0;
		in = sinhx / x;
		ip = (x*coshx-sinhx) / x^2;
		% derivative i0' equal to i1
		din = ip;
	else
		ip = im - (2*n+1)/x*in;
		din = n/(2*n+1)*im + (n+1)/(2*n+1)*ip;
	end
	
	I(n+1,1) = in;
	I(n+1,2) = din;
	
	im = in;
	in = ip;
end


return
end

