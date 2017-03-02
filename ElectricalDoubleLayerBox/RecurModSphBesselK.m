
function [k] = RecurModSphBesselK(nmax,x)

k = recur2(nmax,x);

return
end

function [k] = recur1(nmax,x)
	k = zeros(nmax+1,1);
	
	for n = 0:nmax
		if n == 0
			kn = 1/x;
		elseif n == 1
			kn = (1+x)/x^2;
		else
			kn1 = k(n);
			kn2 = k(n-1);
			kn = (2*n-1)/x*kn1 + kn2;
		end
		k(n+1) = kn;
	end
	
	k = k .* exp(-x);
	
	if 1
		% the sqrt(pi/2) factor
		k = k .* (pi/2);
	end
return
end

function [k] = recur2(nmax,x)
	k = zeros(nmax+1,1);
	
	expnx = exp(-x);
	
	for n = 0:nmax
		if n == 0
			kn = expnx/x;
		elseif n == 1
			kn = expnx/x^2*(x+1);
		else
			kn1 = k(n);
			kn2 = k(n-1);
			kn = (2*n-1)/x*kn1 + kn2;
		end
		k(n+1) = kn;
	end
	
	if 1
		% the sqrt(pi/2) factor
		k = k .* (pi/2);
	end
return
end


