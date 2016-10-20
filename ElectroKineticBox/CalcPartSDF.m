
function [ sdf,rr,rx,ry ] = CalcPartSDF(x0,y0,r0, xs,ys)

rx = xs - x0;
ry = ys - y0;
rr = sqrt(rx.^2 + ry.^2);
sdf = rr - r0;




return
end

