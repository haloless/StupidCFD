
function [S1] = SumS1(x,y,z,lx,ly,lz,kmax)

twopiI = pi * 2.0i;

S1 = 0;

kxs = (-kmax:kmax) ./ lx;
kys = (-kmax:kmax) ./ ly;
kzs = (-kmax:kmax) ./ lz;

[kx,ky] = ndgrid(kxs,kys);
kx2 = kx.^2;
ky2 = ky.^2;
kx_x = kx .* x;
ky_y = ky .* y;

for kz = kzs
	k2 = kx2 + ky2 + kz^2;
	kr = kx_x + ky_y + kz*z;
	ss = exp(twopiI .* kr) ./ k2;
	if kz == 0
		% k2(kmax+1,kmax+1) = 1.0;
		ss(kmax+1,kmax+1) = 0.0;
	end
	S1 = S1 + sum(ss(:));
end


% for kz = kzs
% for ky = kys
% for kx = kxs
	% k2 = kx^2 + ky^2 + kz^2;
	% if k2 > 0
		% kr = kx*x + ky*y + kz*z;
		% ss = exp(twopiI * kr) / k2;
		% S1 = S1 + ss;
	% end
% end
% end
% end

S1 = real(S1);

return
end




