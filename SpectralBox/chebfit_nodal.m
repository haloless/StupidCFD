

function [ ws ] = chebfit_nodal(n,xin,yin)

if isrow(xin)
	xin = xin';
end
if isrow(yin)
	yin = yin';
end

% change to range [-1,1]
xmin = min(xin);
xmax = max(xin);
scale = 2.0 / (xmax-xmin);
xs = (xin-xmin).*scale - 1;

cp = chebpoly1(n, xs);

ws = cp \ yin;


return
end

