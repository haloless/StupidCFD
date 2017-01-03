
clear all;

a = 1;
d = 3;

alpha = acosh(d/a);
c = d / coth(alpha);


xis = linspace(0,alpha, 65);
etas = linspace(0,pi,33);

[xis,etas] = ndgrid(xis,etas);

rs = c .* sin(etas) ./ (cosh(xis)-cos(etas));
zs = c .* sinh(xis) ./ (cosh(xis)-cos(etas));


figure; 
plot(rs(:),zs(:),'.'); 
axis equal;

