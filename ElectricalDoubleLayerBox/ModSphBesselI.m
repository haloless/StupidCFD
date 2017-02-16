function [val] = ModSphBesselI(nu,x)

val = sqrt((pi/2)./x) .* besseli(nu+0.5,x);

return
end

