function [r,theta,phi] = sh_cart2sph(x,y,z)

if nargin == 1
    y = x(2);
    z = x(3);
    x = x(1);
end

r = sqrt(x.^2 + y.^2 + z.^2);

theta = acos(z./r);

phi = atan2(y,x);
if phi < 0
	phi = phi + pi*2;
end

return
end


