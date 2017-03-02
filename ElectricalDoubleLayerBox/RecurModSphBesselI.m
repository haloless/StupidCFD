function [I] = RecurModBesselI(nmax,x)

I = zeros(nmax+1,1);

sinhx = sinh(x);
coshx = cosh(x);

for n = 0:nmax
	if n == 0
		in = sinhx / x;
	elseif n == 1
		in = (x*coshx-sinhx) / x^2;
	else
		in1 = I(n);
		in2 = I(n-1);
		in = in2 - (2*n-1)/x*in1;
	end
	I(n+1) = in;
end


return
end

