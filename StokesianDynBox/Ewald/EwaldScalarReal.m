function [xa,ya,yb,xc,yc] = EwaldScalarReal(xi,r)

xa = 0;
ya = 0;
yb = 0;
xc = 0;
yc = 0;

xir = xi * r;
xir2 = xir^2;
s = r;
s2 = s^2;

erfcxir = erfc(xir);
expxir2 = xi/sqrt(pi) * exp(-xir2);

ya = (0.75 + 0.5 / s2) / s * erfcxir + ((1.0 + xir2 * (14.0 + 4.0 * xir2 *  (- 5.0 + xir2))) / s2 - 4.5 + 3.0 * xir2) * expxir2;
a2 = (0.75 - 1.5 / s2) / s * erfcxir + ((- 3.0 + xir2 * (- 2.0 + 4.0 * xir2 * (4.0 - xir2))) / s2 + 1.5 - 3.0 * xir2) * expxir2;
xa = a2 + ya;

yb = - 0.75 / s2 * erfcxir - 1.5 * (+ 1.0 + xir2 * (- 6.0 + xir2 * (+ 2.0))) / s * expxir2;

yc = - 3.0 / 8.0 / s2 / s * erfcxir - 0.75 * (+ 1.0 + xir2 * (+ 14.0 + xir2 * (-20.0 + xir2 * ( + 4.0)))) / s2 * expxir2;
c2 = 9.0 / 8.0 / s2 / s * erfcxir - 0.75 * (- 3.0 + xir2 * (- 2.0 + xir2 * (+ 16.0 + xir2 * (- 4.0)))) / s2 * expxir2;
xc = c2 + yc;




return
end



