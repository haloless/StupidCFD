
function [k] = RecurModSphBesselK2(nmax,x)

k = recur2(nmax,x);

return
end

function [k] = recur2(nmax,x)
	k = zeros(nmax+1,2);
	
	expnx = exp(-x);
	
	% kn = expnx / x;
	% kp = expnx / x^2 * (x+1);
	
	for n = 0:nmax
		if n == 0
			% initial value for n=0,1
			km = 0;
			kn = expnx / x;
			kp = expnx / x^2 * (x+1);
			% derivative n=0 is -k1
			dkn = -kp;
		else
			kp = km + (2*n+1)/x*kn;
			dkn = -n/(2*n+1)*km - (n+1)/(2*n+1)*kp;
		end
		
		k(n+1,1) = kn;
		k(n+1,2) = dkn;
		
		km = kn;
		kn = kp;
	end
	
	if 1
		% the sqrt(pi/2) factor
		k = k .* (pi/2);
	end
return
end


