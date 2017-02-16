function [p] = LegendrePoly(n,x)

% associate Legendre function
val = legendre(n,x);

p = val(1,:);

p = p';

return
end
