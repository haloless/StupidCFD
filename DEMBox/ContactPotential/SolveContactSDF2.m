function [pos,ok] = SolveContactSDF2(sdf1,sdf2, xguess)

% length scale
lref = (sdf1.xmax-sdf1.xmin) / 10;

% scale SDF value at O(1)
% sdfscale = 1.0 / lref;
sdfscale = 1.0;

sdfmin1 = min(sdf1.phig(:));
sdfmin2 = min(sdf2.phig(:));
sdfmin1 = -min(sdfmin1,sdfmin2);
sdfmin2 = sdfmin1;

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
		% use mu^2 = (1-lambda)/(1+lambda)
		rhs = -[ g1+lambda^2*g2; sdfscale*(p1-p2) ];
		
		amat = h1 + lambda^2*h2;
		bvec = 2*lambda*g2;
		cvec = sdfscale*(g1-g2);
		mat = [ amat, bvec; cvec', 0 ];
	end
	if 1
		% use mu^2 = (1-lambda)/(1+lambda)
		
		P1 = p1 * (p1+2*sdfmin1);
		G1 = 2*(p1+sdfmin1)*g1;
		H1 = 2*g1*g1' + 2*(p1+sdfmin1)*h1;
		
		P2 = p2 * (p2+2*sdfmin2);
		G2 = 2*(p2+sdfmin2)*g2;
		H2 = 2*g2*g2' + 2*(p2+sdfmin2)*h2;
		
		rhs = -[ G1+lambda^2*G2; sdfscale*(P1-P2) ];
		
		amat = H1 + lambda^2 * H2;
		bvec = 2*lambda * G2;
		cvec = sdfscale*(G1-G2);
		mat = [ amat, bvec; cvec', 0 ];
	end
	
	
	
	rnorm = norm(rhs);
	disp(['iter=',int2str(iter), '; res=',num2str(rnorm), ...
    '; rhs=',num2str(rhs'), '; x=',num2str(x')]);
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

% contact flag
ok = 0;
if p1<0 && p2<0
	ok = 1;
end

return
end







