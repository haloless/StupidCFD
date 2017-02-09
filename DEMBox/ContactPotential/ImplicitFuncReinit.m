function [phi] = ImplicitFuncReinit(phi0,nx,ny,dx,dy, bandwidth)
% Level-set reinitialization.
%

cfl = 0.2;
dt = cfl * dx;
niter = ceil(bandwidth/dt);

% initial sign
sgn0 = CalcSign(phi0,nx,ny,dx,dy);

phi = phi0;

for iter = 1:niter
	
	% 2nd-order RK
	
	phi1 = Reinit(phi0,sgn0,phi,nx,ny,dx,dy,dt);
	
	phi2 = Reinit(phi0,sgn0,phi1,nx,ny,dx,dy,dt);
	
	phin = (phi + phi2) .* 0.5;
	
	phi = phin;
	
	if mod(iter,50) == 0
		fprintf('iter=%d in %d\n', iter,niter);
	end
end



return
end

function [sgn] = CalcSign(phi,nx,ny,dx,dy)

dh = dx;

% smeared sign
sgn = phi ./ sqrt(phi.^2 + dh^2);

if 0
	for i = 1:nx
	for j = 1:ny
		if abs(phi(i,j)) < dh
			sgn(i,j) = 0;
		end
	end
	end
end

return
end

function [phinew] = Reinit(phi0,sgn0,phi,nx,ny,dx,dy,dt)

phinew = phi;

for j = 1:ny
for i = 1:nx
	
	dmx = 0;
	if i > 1
		dmx = (phi(i,j)-phi(i-1,j)) / dx;
	end
	
	dpx = 0;
	if i < nx
		dpx = (phi(i+1,j)-phi(i,j)) / dx;
	end
	
	dmy = 0;
	if j > 1
		dmy = (phi(i,j)-phi(i,j-1)) / dy;
	end
	
	dpy = 0;
	if j < ny
		dpy = (phi(i,j+1)-phi(i,j)) / dy;
	end
	
	s = sgn0(i,j);
	if s > 0
		fx = max(abs(max(dmx,0)), abs(min(dpx,0)));
		fy = max(abs(max(dmy,0)), abs(min(dpy,0)));
	else
		fx = max(abs(max(dpx,0)), abs(min(dmx,0)));
		fy = max(abs(max(dpy,0)), abs(min(dmy,0)));
	end
	
	% |grad phi|
	gphi = sqrt(fx^2+fy^2);
	
	phinew(i,j) = phi(i,j) + dt*s*(1-gphi);
end
end

return
end


