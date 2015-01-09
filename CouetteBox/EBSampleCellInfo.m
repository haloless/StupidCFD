
function [ frac uint vint ] = EBSampleCellInfo(xcen,ycen,dx,dy)

% Couette flow, inner cylinder
global Omega0 R0 R1

nsub = 5;

if (1)

[xx,yy] = ndgrid( ...
linspace(xcen-0.5*dx+0.5*dx/nsub, xcen+0.5*dx-0.5*dx/nsub, nsub), ...
linspace(ycen-0.5*dy+0.5*dy/nsub, ycen+0.5*dy-0.5*dy/nsub, nsub));

dist = ProbDistFunc(xx,yy);
mask = (dist <= 0);

%
frac = double(sum(mask(:))) / nsub^2;

%
rr = sqrt(xx.^2 + yy.^2);
uu = -Omega0 .* yy .* (rr < (R0+R1)/2);
vv =  Omega0 .* xx .* (rr < (R0+R1)/2);
uint = uu .* double(mask);
uint = sum(uint(:)) / nsub^2;
vint = vv .* double(mask);
vint = sum(vint(:)) / nsub^2;

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


