
function [xs,scale,xmin,xmax] = chebrange(xin)

xmin = min(xin);
xmax = max(xin);
scale = 2.0 / (xmax-xmin);
xs = (xin-xmin).*scale - 1;

return
end

