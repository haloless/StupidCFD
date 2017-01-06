function [pos] = SolveContactPotential(shapea,shapeb)

if shapea.type~=2 || shapeb.type~=2
	error('Must use SuperEllipse');
end

disp(mfilename());

tmpa = CalcShapeInfo(shapea);
tmpb = CalcShapeInfo(shapeb);

xguess = (tmpa.xc+tmpb.xc)/2;
if max(tmpa.p)>2
	xguess = SolveIntersectPotential(shapea,tmpa.xc,tmpb.xc-tmpa.xc);
end
x = xguess;
lambda = 0;

maxiter = 25;
tolrel = 1.0e-6;
tolabs = 1.0e-9;

conv = 0;

% initial residual
ga = CalcG(tmpa,x);
gb = CalcG(tmpb,x);
ha = CalcH(tmpa,x);
hb = CalcH(tmpb,x);
resid = [gb+lambda*ga; ShapePotential(shapea,x(1),x(2))];
rnorm0 = norm(resid);
maxincr = min(tmpa.a);

disp(['|res0|=',num2str(rnorm0)]);

for iter = 1:maxiter
	
	mat = [hb+lambda*ha, ga; ga', 0];
	
	rhs = -resid;
	
	% sol = mat \ rhs;
	invmat = SafeInvert(mat);
	sol = invmat * rhs;
	
	omega = 1.0;
	xincr = sol(1:2);
	if norm(xincr) > maxincr
		omega = maxincr / norm(xincr);
	end
	
	x = x + xincr * omega;
	lambda = lambda + sol(3) * omega;
	
	ga = CalcG(tmpa,x);
	gb = CalcG(tmpb,x);
	ha = CalcH(tmpa,x);
	hb = CalcH(tmpb,x);
	resid = [gb+lambda*ga; ShapePotential(shapea,x(1),x(2))];
	
	rnorm = norm(resid);
	disp(['iter=',int2str(iter), '; |res|=',num2str(rnorm)]);
	
	if rnorm < rnorm0*tolrel+tolabs
		conv = 1;
		break;
	end
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









