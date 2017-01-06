function [pos] = SolveContactWall(shape,wall)

if shape.type~=2
	error('Must use SuperEllipse');
end

disp(mfilename());

tmp = CalcShapeInfo(shape);

% xguess = (tmp.xc+wall.xwall)/2;
xguess = tmp.xc - mean(tmp.a)*wall.nwall;
if max(tmp.p) > 2
	xguess = SolveIntersectPotential(shape, tmp.xc,-wall.nwall);
end
x = xguess;
lambda = 0.0;

maxiter = 25;
tolrel = 1.0e-6;
tolabs = 1.0e-9;
maxincr = 0.5*(min(tmp.a));

conv = 0;

% initial residual
ga = CalcG(tmp,x);
ha = CalcH(tmp,x);
gw = wall.nwall;
hw = zeros(2,2);
% resid = [gw+lambda*ga; ShapePotential(shape,x(1),x(2))];
% resid = [ga+lambda*gw; WallPotential(wall,x(1),x(2))];
resid = [ gw(1)*ga(2)-gw(2)*ga(1); ShapePotential(shape,x(1),x(2)) ];
rnorm0 = norm(resid);
disp(['|res0|=',num2str(rnorm0)]);
% if rnorm0 < tolabs
	% pos = x;
	% return;
% end


for iter = 1:maxiter
	
	% mat = [hw+lambda*ha, ga; ga', 0];
	% mat = [ha+lambda*hw, gw; gw', 0];
	mat = [ gw(1)*ha(1,2)-gw(2)*ha(1,1), gw(1)*ha(2,2)-gw(2)*ha(1,2); ga(1),ga(2) ];
	
	rhs = -resid;
	
	% sol = mat \ rhs;
	minv = SafeInvert(mat);
	sol = minv * rhs;
	
	omega = 1.0;
	xincr = sol(1:2);
	if norm(xincr) > maxincr
		omega = maxincr / norm(xincr);
	end
	
	x = x + xincr*omega;
	% lambda = lambda + sol(3);
	
	ga = CalcG(tmp,x);
	ha = CalcH(tmp,x);
	% resid = [gw+lambda*ga; ShapePotential(shape,x(1),x(2))];
	% resid = [ga+lambda*gw; WallPotential(wall,x(1),x(2))];
	resid = [ gw(1)*ga(2)-gw(2)*ga(1); ShapePotential(shape,x(1),x(2)) ];
	
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



