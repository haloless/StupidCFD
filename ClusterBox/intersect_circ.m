
function [isec,dsec] = intersect_circ(a,x,y,u,v,len, np,xs,ys)

isec = 0;
dsec = 0;

asmall = a * 1.0e-5;

I = 1:np;
xrel = x - xs(I);
yrel = y - ys(I);

bs = u.*xrel + v.*yrel;
cs = xrel.^2 + yrel.^2 - a^2;

disc = bs.^2 - cs;
ids = find(disc >= 0);
nprob = numel(ids);
if nprob == 0
	return
end

dist1 = -bs(ids) - sqrt(disc(ids));
dist2 = -bs(ids) + sqrt(disc(ids));

ok = ones(nprob,1);
dist = zeros(nprob,1);
for i = 1:nprob
	if abs(dist1(i)) < abs(dist2(i))
		dist(i) = dist1(i);
	else
		dist(i) = dist2(i);
	end
	
	if len >= 0
		if dist(i)>=-asmall && abs(dist(i))<=len
		% if abs(dist(i))<=len
			ok(i) = 1;
		else
			ok(i) = 0;
		end
	end
end

isec = 0;
dsec = 99999;
for i = 1:nprob
	if ok(i) == 1
		if abs(dist(i)) < abs(dsec)
			isec = ids(i);
			dsec = dist(i);
		end
	end
end


return
end

