function [Bnm] = FuncB(n,m,kR)

numax = min(n,m);

Bnm = 0;
for nu = 0:numax
	Bnm = Bnm + FuncA(nu,n,m) * ModSphBesselK(n+m-2*nu,kR);
end

return
end

