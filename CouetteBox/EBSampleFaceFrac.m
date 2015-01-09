
function [ frac ] = EBSampleFaceFrac(dir,xcen,ycen,dx,dy)

nsub = 25;

if dir == 0
    xx = zeros(nsub,1) + xcen;
    yy = linspace(ycen-0.5*dy+0.5*dy/nsub, ycen+0.5*dy-0.5*dy/nsub, nsub)';
else
    xx = linspace(xcen-0.5*dx+0.5*dx/nsub, xcen+0.5*dx-0.5*dx/nsub, nsub)';
    yy = zeros(nsub,1) + ycen;
end

dist = ProbDistFunc(xx,yy);
frac = double(sum(dist(:)<=0)) / nsub;

return
end


