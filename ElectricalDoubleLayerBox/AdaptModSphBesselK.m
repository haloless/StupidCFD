function [khat] = AdaptModSphBesselK(n,x)

coef = exp(x) .* (x.^(n+1)) ./ double_factorial(2*n-1);

kn = ModSphBesselK(n,x);
khat = 2/pi .* coef .* kn;

return
end

