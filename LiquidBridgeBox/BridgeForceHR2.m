function [F] = BridgeForceHR2(R1,R2,H,theta1,theta2,V,sigma)

if ~exist('sigma','var')
	sigma = 1.0;
end

if R2 < 0
	R2 = 0;
end

if 0
	% initial guess with straight cylindrical profile
	[~,alpha1] = SolveStraightBridge(R1,R2,H,V);
else
	alpha1 = estim_alpha(R1,R2,H,V);
end

% this is the volume function V = V(alpha1)
% it is to solve the embracing angle alpha1 with this nonlinear relationship
vfun = @(a1) calc_vol(R1,R2,H,theta1,theta2,a1) - V;

disp(mfilename);
if 0
	% use Matlab FSOLVE
	alpha1 = fsolve(vfun,alpha1);
else
	% use our own Newton's method
	% tolerance is based on the input volume
	tol = V * 1.0e-4;
	alpha1 = SolveFunc(vfun,alpha1,tol);
end

alpha2 = calc_alpha2(R1,R2,H,theta1,theta2, alpha1);

[rho1,rho2] = calc_rho(R1,R2,H,theta1,theta2, alpha1,alpha2);

% check neck
xneck = R1*cos(alpha1) + rho1*cos(alpha1+theta1);
necklb = R1*cos(alpha1);
neckub = R1+H+R2 - R2*cos(alpha2);
if (xneck<necklb || xneck>neckub)
	b1 = R1*sin(alpha1);
	if R2 > 0
		b2 = R2*sin(alpha2);
	else
		b2 = b1 + rho1*sin(alpha1+theta1) - rho1*sin(theta2);
	end
	% rho2 = 0.5*(b1+b2);
	
	rho2a = b1 / sin(alpha1+theta1);
	rho2b = b2 / sin(alpha2+theta2);
	rho2 = 0.5*(rho2a+rho2b);
	
	pres = sigma * (1/rho2-1/rho1);
	F1 = AxisymEvalForce(alpha1,theta1,b1,pres,sigma);	
	F2 = AxisymEvalForce(alpha2,theta2,b2,pres,sigma);
	
	F = (F1+F2)/2;
else
	F = pi*sigma * rho2*(1+rho2/rho1);
end


return
end

function [alpha2] = calc_alpha2(R1,R2,H,theta1,theta2, alpha1)
	if R2 > 0
		A = tan(alpha1/2);
		C = tan((theta1-theta2)/2);
		By = (R1*A + H*A/2 + H*C/2);
		Bx = (R2+H/2-R1*A*C-R2*A*C-H*A*C/2);
		alpha2 = 2 * atan2(By,Bx);
	else
		% wall
		alpha2 = 0;
	end
	return
end

function [rho1,rho2] = calc_rho(R1,R2,H,theta1,theta2, alpha1,alpha2)
	at1 = alpha1 + theta1;
	at2 = alpha2 + theta2;
	
	rho1 = (R1*(1-cos(alpha1)) + R2*(1-cos(alpha2)) + H) / (cos(at1)+cos(at2));
	if 1
		% limit outer radius
		rho1max = R1 * 1000;
		if abs(rho1) > rho1max
			rho1 = rho1max * sign(rho1);
		end
	end
	
	
	rho2 = R1*sin(alpha1) - rho1*(1-sin(at1));
	
	return
end

function [vol] = calc_vol(R1,R2,H,theta1,theta2, alpha1)
	
	alpha2 = calc_alpha2(R1,R2,H,theta1,theta2, alpha1);
	
	at1 = alpha1 + theta1;
	at2 = alpha2 + theta2;
	ca1 = cos(alpha1);
	ca2 = cos(alpha2);
	cat1 = cos(at1);
	cat2 = cos(at2);
	sat1 = sin(at1);
	sat2 = sin(at2);
	
	[rho1,rho2] = calc_rho(R1,R2,H,theta1,theta2, alpha1,alpha2);
	
	vrot = ((rho1+rho2)^2*rho1+rho1^3) * (cat1+cat2);
	vrot = vrot + (rho1+rho2)*rho1^2*(at1+at2-pi);
	vrot = vrot - (rho1+rho2)*rho1^2*(sat1*cat1+sat2*cat2);
	vrot = vrot - 1.0/3.0*rho1^3*(cat1^3+cat2^3);
	vrot = vrot * pi;
	
	vcap1 = pi/3 * R1^3 * (2-3*ca1+ca1^3);
	vcap2 = pi/3 * R2^3 * (2-3*ca2+ca2^3);
	
	vol = vrot - vcap1 - vcap2;
	return
end

function [alpha] = estim_alpha(R1,R2,H,V)
	% if R2 > 0
		% R = DerjaguinRadius(R1,R2);
		% h = -H/2 + 0.5*sqrt(H^2+2*V/pi/R);
	% else
		% R = R1;
		% h = -H + sqrt(H^2+V/pi/R);
	% end
	% alpha = sqrt(2*h/R);
	
	if R2 > 0
		coef = 1 + R1/R2;
	else
		coef = 1;
	end
	
	Hstar = H / R1;
	Vstar = V / (pi*R1^3);
	a2 = (-Hstar + sqrt(Hstar^2+coef*Vstar)) / coef * 2;
	alpha = sqrt(a2);
	
	return
end




