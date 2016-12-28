
function [u] = SumU(x,y,z,lx,ly,lz,kmax)

twopiI = -pi * 2.0i;

u = 0;

% kxs = (-kmax:kmax) ./ lx;
kxs = (0:kmax) ./ lx;
kys = (-kmax:kmax) ./ ly;
kzs = (-kmax:kmax) ./ lz;

[ky,kz] = ndgrid(kys,kzs);
% ky = reshape(ky,[],1);
% kz = reshape(kz,[],1);

ky2 = ky.^2;
kz2 = kz.^2;
ky2kz2 = ky2 + kz2;

kz_z = kz .* z;
ky_y = ky .* y;

eyz = exp(twopiI .* (ky_y+kz_z));
eyz = real(eyz);

cyz = -ky2kz2 .* eyz;

for kx = kxs
	k2 = kx.^2 + ky2kz2;
	k4 = k2.^2;
	
	% uu = cyz ./ k4 .* exp(twopiI*kx*x);
	
	cx = 2;
	if kx == 0
		cx = 1;
	end
	uu = cyz ./ k4 .* (cx*cos(2*pi*kx*x));

	if kx == 0
		uu(kmax+1,kmax+1) = 0;
	end
	
	u = u + sum(uu(:));
end

% u = real(u);

fx = 1.0;
mu = 1.0/(6*pi);
u = u * (-fx) / (4*pi^2*mu*lx*ly*lz);

return
end

