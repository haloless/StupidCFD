function [pos] = SolveIntersectPotential(shape, xbase,evec)
%% x = xbase + evec * t
%% find t

if shape.type~=2
	error('Must use SuperEllipse');
end

disp(mfilename());

tmp = CalcShapeInfo(shape);

% use shape size scale as first guess
tguess = mean(tmp.a);
t = tguess;
x = xbase + t * evec;

maxiter = 25;
tolrel = 1.0e-6;
tolabs = 1.0e-12;

conv = 0;

% initial residual
ga = CalcG(tmp,x);
ha = CalcH(tmp,x);
resid = [ ShapePotential(shape,x(1),x(2)) ];
rnorm0 = norm(resid);
disp(['|res0|=',num2str(rnorm0)]);

for iter = 1:maxiter
	
	mat = [ dot(ga,evec) ];
	
	rhs = -resid;
	
	sol = mat \ rhs;
	
	t = t + sol(1);
	x = xbase + t*evec;
	
	ga = CalcG(tmp,x);
	ha = CalcH(tmp,x);
	resid = [ ShapePotential(shape,x(1),x(2)) ];
	
	rnorm = norm(resid);
	disp(['iter=',int2str(iter), '; |res|=',num2str(rnorm)]);
	
	if rnorm < rnorm0*tolrel+tolabs
		conv = 1;
		break;
	end
end

if ~conv
	% error('Failed');
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


