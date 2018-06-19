function [p] = LegendrePoly(n,x)

% associate Legendre function
val = legendre(n,x);

% the 0th-order assoc Legendre is Legendre Poly.
p = val(1,:);

p = p';

return
end
