
function [ frac ] = CalcPartFrac(x0,y0,r0, xs,ys,dx,dy)

[sdf,rr,rx,ry] = CalcPartSDF(x0,y0,r0, xs,ys);

frac = zeros(size(xs));
cutoff = max(dx,dy) * 2;
frac(sdf>=cutoff) = 0;
frac(sdf<=-cutoff) = 1;

range = find(sdf>-cutoff & sdf<cutoff)';
for ind = range
    nvx = rx(ind) / rr(ind);
    nvy = ry(ind) / rr(ind);
    dist = -sdf(ind);
    
    frac(ind) = CubeChopGetVolume2D(nvx,nvy,dist, dx,dy);
end



return
end

