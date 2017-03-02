function [ihat] = AdaptModSphBesselI(n,x)

coef = double_factorial(2*n+1) ./ (x.^n);

in = ModSphBesselI(n,x);

ihat = coef .* in;

return
end

