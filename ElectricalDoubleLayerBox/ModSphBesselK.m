function [val] = ModSphBesselK(nu,x)

val = sqrt((2/pi)./x) .* besselk(nu+0.5,x);

return
end

