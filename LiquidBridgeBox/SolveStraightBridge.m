%% 
%% Assume the bridge is a straight cylinder connecting two spheres.
%% Calculate the cylinder radius for given volume.
%%
function [r,alpha1,alpha2] = SolveStraightBridge(R1,R2,H,V)

Vinv = 1/V;

% normalize by input volume for numerical reason
vfunc = @(a) StraightVolume(R1,R2,H,a).*Vinv - 1;

% solve cylinder radius
disp(mfilename);
rguess = R1*0.1;
if 1
	% use Matlab FSOLVE
	options = optimoptions('fsolve','TolFun',1.0e-3,'TolX',1.0e-8);
	r = fsolve(vfunc,rguess,options);
else
	% use our own Newton's method
	% tolerance based on input volume
	r = SolveFunc(vfunc,rguess,1.0e-3);
end

% Vsol = StraightVolume(R1,R2,H,r)

% calculate two embracing angles
alpha1 = asin(r/R1);
if R2 > 0
	alpha2 = asin(r/R2);
else
	alpha2 = 0;
end

return
end

function [vol] = StraightVolume(R1,R2,H,r)
	
	alpha1 = asin(r/R1);
	h1 = R1 * (1-cos(alpha1));
	vol1 = SphereCapVolume(r,h1);
	
	if R2 > 0
		alpha2 = asin(r/R2);
		h2 = R2 * (1-cos(alpha2));
		vol2 = SphereCapVolume(r,h2);
	else
		alpha2 = 0;
		h2 = 0;
		vol2 = 0;
	end
	
	
	
	% cylinder volume
	vol = pi*r^2*(H+h1+h2);
	% minus two sphere cap
	vol = vol - vol1 - vol2;
	
	return
end


