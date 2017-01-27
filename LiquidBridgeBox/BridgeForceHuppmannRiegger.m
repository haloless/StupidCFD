%% 
%% The bridge force equation of (Huppmann & Riegger, 1975).
%%
function [F] = BridgeForceHuppmannRiegger(R,H,theta,V,sigma)

if ~exist('sigma','var')
	sigma = 1.0;
end


%
vfun = @(alpha) vol_func(R,H,theta,alpha)-V;

% alpha0 = pi/6;
[~,alpha0] = SolveStraightBridge(R,R,H,V);

disp(mfilename);
alpha = fsolve(vfun,alpha0);

%
[rho1,rho2] = rho_func(R,H,theta,alpha);

F = 2*pi*sigma*R*sin(alpha)*sin(alpha+theta) + pi*sigma*(R*sin(alpha))^2*(1/rho2-1/rho1);

return
end

%% The inner and outer curvature radius.
function [rho1,rho2] = rho_func(R,H,theta,alpha)
	sina = sin(alpha);
	cosa = cos(alpha);
	sinat = sin(alpha+theta);
	cosat = cos(alpha+theta);
	
	rho1 = R*sina - (R*(1-cosa)+H/2)*(1-sinat)/cosat;
	rho2 = (R*(1-cosa)+H/2) / cosat;
	return
end

function [vol] = vol_func(R,H,theta,alpha)
	cosat = cos(alpha+theta);
	
	[rho1,rho2] = rho_func(R,H,theta,alpha);
	
	vol = 2*pi*(cosat-(pi/2-alpha-theta))*(rho2^3+rho1*rho2^2) + pi*rho1^2*rho2*cosat;
return
end


