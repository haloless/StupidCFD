function [acoef] = SphHarmReexpandApply(nmax,Rmat,Smat,bcoef)


xtmp = zeros(size(bcoef));
for n = 0:nmax
	for m = -n:n
		inm = sh_sub2ind(n,m);
		xnm = 0;
		for s = -n:n
			ins = sh_sub2ind(n,s);
			xnm = xnm + Rmat(inm,ins) * bcoef(ins);
		end
		xtmp(inm) = xnm;
	end
end

ytmp = zeros(size(bcoef));
for n = 0:nmax
	for m = -n:n
		inm = sh_sub2ind(n,m);
		ynm = 0;
		for l = abs(m):nmax
			ilm = sh_sub2ind(l,m);
			ynm = ynm + Smat(inm,ilm) * xtmp(ilm);
		end
		ytmp(inm) = ynm;
	end
end

ztmp = zeros(size(bcoef));
for n = 0:nmax
	for m = -n:n
		inm = sh_sub2ind(n,m);
		znm = 0;
		for s = -n:n
			ins = sh_sub2ind(n,s);
			znm = znm + conj(Rmat(ins,inm)) * ytmp(ins);
		end
		ztmp(inm) = znm;
	end
end

acoef = ztmp;

return
end

