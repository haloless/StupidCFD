
function [ frac ] = EBSampleCellFrac(xcen,ycen,dx,dy)

% nsub = 5;
nsub = 8;

if (1)
[xx,yy] = ndgrid( ...
linspace(xcen-0.5*dx+0.5*dx/nsub, xcen+0.5*dx-0.5*dx/nsub, nsub), ...
linspace(ycen-0.5*dy+0.5*dy/nsub, ycen+0.5*dy-0.5*dy/nsub, nsub));

dist = ProbDistFunc(xx,yy);
frac = double(sum(dist(:)<=0)) / nsub^2;

else
xlo = xcen - dx*0.5;
ylo = ycen - dy*0.5;


frac = 0.0;

for isub = 1:nsub
for jsub = 1:nsub
    xx = xlo + (isub-0.5)*(dx/nsub);
    yy = ylo + (jsub-0.5)*(dy/nsub);
    if ProbDistFunc(xx, yy) <= 0
        frac = frac + 1.0;
    end
end
end

frac = frac / nsub^2;

end

return
end


