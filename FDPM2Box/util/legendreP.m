function [p] = legendreP(n, x)

ps = legendre(n, x);

p = ps(1,:);

p = reshape(p, size(x));

return
end


