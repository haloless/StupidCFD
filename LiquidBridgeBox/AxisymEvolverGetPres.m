%% Pressure is the Lagrange multiplier for the volume constraint.
%% Sometimes, Matlab FMINCON may return negative value.
%% This is possible for nonlinear equality constraint.
%% Therefore in this function we check the correct pressure.
%% The current way is to compute the curvature.
function [pres] = AxisymEvolverGetPres(bridge, np,rp,xp, lambda)


R1 = bridge.R1;
R2 = bridge.R2;
H = bridge.H;
theta1 = bridge.theta1;
theta2 = bridge.theta2;
sigma = bridge.sigma;

% Lagrange multiplier
pres = lambda.eqnonlin(3);

if 0
	%
	r1 = rp(1);
	r2 = rp(np);
	alpha1 = AxisymEmbraceAngle(R1,r1);
	alpha2 = AxisymEmbraceAngle(R2,r2);

	% with positive pressure
	[F1a] = AxisymEvalForce(alpha1,theta1,r1,pres,sigma);
	[F2a] = AxisymEvalForce(alpha2,theta2,r2,pres,sigma);

	% with negative pressure
	[F1b] = AxisymEvalForce(alpha1,theta1,r1,-pres,sigma);
	[F2b] = AxisymEvalForce(alpha2,theta2,r2,-pres,sigma);

	if abs(F1a-F2a) <= abs(F1b-F2b)
		pres = pres;
	else
		pres = -pres;
	end
end

if 1
	% evaluate curvature near middle part of the bridge
	% although the curvature will not be very accurate
	% its sign can be used to fix the Laplace pressure
	ielem = round(np/4);
	[curv,ko,ki] = EstimCurv(bridge, np,rp,xp, ielem);
	pres = abs(pres) * sign(curv);
end


return
end

function [curv,ko,ki] = EstimCurv(bridge, np,rp,xp, ielem)
	
	if (ielem==1 || ielem==np-1)
		error('Estimation based on center difference, not applicable to two end elements!');
	end
	
	xs = xp(ielem-1:ielem+2);
	rs = rp(ielem-1:ielem+2);
	
	dx = xs(2:4) - xs(1:3);
	dr = rs(2:4) - rs(1:3);
	dl = sqrt(dx.^2 + dr.^2);
	
	tx = dx ./ dl;
	tr = dr ./ dl;
	tvec = [tx, tr]';
	
	% tangent vector at two ends of center element
	tl = tvec(:,1) + tvec(:,2);
	tl = tl ./ norm(tl);
	tr = tvec(:,2) + tvec(:,3);
	tr = tr ./ norm(tr);
	
	% tangent and normal vector at center element
	tc = tvec(:,2);
	nc = [-tc(2); tc(1)];
	dlc = dl(2);
	
	% finite difference, dt/dl
	dtds = (tr-tl)./dlc;
	% signed curvature of outer surface
	% note it has the same sign of 2nd derivative of the surface curve
	% so ko>0 for concave surface, and ko<0 for convex
	ko =  dot(dtds,nc);
	
	angc = atan2(tc(2),tc(1));
	rc = 0.5*(rs(2)+rs(3));
	
	% curvature of inner side
	% this should be always > 0
	ki = rc / cos(angc);
	ki = 1/ki;
	
	curv = ki - ko;
	
	return
end


