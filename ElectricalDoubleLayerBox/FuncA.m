function [A] = FuncA(nu,n,m)

a1 = gamma(n-nu+1/2) .* gamma(m-nu+1/2) .* gamma(nu+1/2) ./ gamma(m+n-nu+3/2);

a2 = factorial(n+m-nu) ./ factorial(n-nu) ./ factorial(m-nu) ./ factorial(nu);

A = (n+m-2*nu+1/2)/pi .* a1 .* a2;

return
end
