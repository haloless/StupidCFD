
function [gl] = FuncGReg(rr,xx,alpha,xi)

if rr > 0
	er = xi * rr;
	ar = alpha * rr;
	
	c1 = (erf(er) - erf(ar)) / rr;
	c2 = 2.0/sqrt(pi) * (xi*exp(-er*er) - alpha*exp(-ar*ar));
	
	ee = xx ./ rr;
	ee = ee * ee';
	
	ee1 = ee + eye(3);
	ee2 = -ee + eye(3);
	
	gl = ee1*c1 + ee2*c2;
else
	gl = eye(3) .* (4*(xi-alpha)/sqrt(pi));
end


return
end

