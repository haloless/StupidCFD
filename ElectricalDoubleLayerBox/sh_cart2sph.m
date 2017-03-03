function [r,theta,phi] = sh_cart2sph(x,y,z)

r = sqrt(x.^2 + y.^2 + z.^2);

theta = acos(z./r);

phi = atan2(y,x);
if phi < 0
	phi = phi + pi*2;
end

return
end


