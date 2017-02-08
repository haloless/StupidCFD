function [pos] = SolveContactSDF(sdf1,sdf2, xguess)

% length scale
lref = (sdf1.xmax-sdf1.xmin) / 10;

% scale SDF value at O(1)
sdfscale = 1.0 / lref;

disp(mfilename());

x = xguess;
lambda = 1.0;

maxiter = 25;
tolrel = 1.0e-3;
tolabs = 1.0e-3;

conv = 0;


for iter = 1:maxiter
	
	[p1,g1,h1] = SDFPotential(sdf1, x(1),x(2));
	[p2,g2,h2] = SDFPotential(sdf2, x(1),x(2));
	
	if 0
		% use lambda 
		rhs = -[ (1+lambda)*g1 + (1-lambda)*g2; sdfscale*(p1-p2) ];
		
		amat = (1+lambda)*h1 + (1-lambda)*h2;
		bvec = sdfscale*(g1-g2);
		
		mat = [ amat, bvec; bvec', 0 ];
	end
	if 1
		% use mu^2 = (1-lambda)/(1+lambda)
		rhs = -[ g1+lambda^2*g2; sdfscale*(p1-p2) ];
		
		amat = h1 + lambda^2*h2;
		bvec = 2*lambda*g2;
		cvec = sdfscale*(g1-g2);
		mat = [ amat, bvec; cvec', 0 ];
	end
	if 0
		% use mu = (1-lambda)/(1+lambda)
		rhs = -[ g1+lambda*g2; sdfscale*(p1-p2) ];
		
		amat = h1 + lambda*h2;
		bvec = 2*g2;
		cvec = sdfscale*(g1-g2);
		mat = [ amat, bvec; cvec', 0 ];
	end
	
	
	
	rnorm = norm(rhs);
	disp(['iter=',int2str(iter), '; res=',num2str(rnorm), '; rhs=',num2str(rhs')]);
	if iter == 1
		rnorm0 = rnorm;
	end
	if rnorm <= rnorm0*tolrel+tolabs
		conv = 1;
		break;
	end
	
	sol = mat \ rhs;
	% invmat = SafeInvert(mat);
	% sol = invmat * rhs;
	
	omega = 1.0;
	xincr = sol(1:2);
	% if norm(xincr) > maxincr
		% omega = maxincr / norm(xincr);
	% end
	
	x = x + xincr * omega;
	lambda = lambda + sol(3) * omega;
end

if ~conv
	error('Failed');
end

pos = x;

return
end

function [tmp] = CalcShapeInfo(shape)
	tmp = struct();
	
	tmp.a = [shape.a; shape.b];
	tmp.p = [shape.p; shape.q];
	tmp.pa = tmp.p ./ tmp.a;
	tmp.p1a = (tmp.p-1) ./ tmp.a;
	
	tmp.R = RotationMatrix(-shape.rotang);
	tmp.xc = [shape.xc; shape.yc];
	
	return
end

% gradient
function [g] = CalcG(tmp,x)
	%
	x0 = tmp.R * (x-tmp.xc);
	
	y0 = (x0./tmp.a).^(tmp.p-1);
	
	gvec = tmp.pa .* y0;
	
	g = tmp.R' * gvec;
	
	return
end

% hessian
function [h] = CalcH(tmp,x)
	%
	x0 = tmp.R * (x-tmp.xc);
	
	y0 = (x0./tmp.a).^(tmp.p-2);
	
	hvec = tmp.pa .* tmp.p1a .* y0;
	
	h = tmp.R' * diag(hvec) * tmp.R;
	
return
end









