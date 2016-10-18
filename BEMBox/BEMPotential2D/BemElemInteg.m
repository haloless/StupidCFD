
function [ aa,bb ] = BemElemInteg(i,j)
% A(i,j) and b(i)

BemMeshGlobals;

al = dlen(j);

j1 = node(1,j);
j2 = node(2,j);

x11 = y(1,j1) - x(1,i);
x21 = y(2,j1) - x(2,i);
x12 = y(1,j2) - x(1,i);
x22 = y(2,j2) - x(2,i);

r1 = sqrt(x11^2 + x21^2);
r2 = sqrt(x12^2 + x22^2);
d = x11*dnorm(1,j) + x21*dnorm(2,j);
t1 = -x11*dnorm(2,j) + x21*dnorm(1,j);
t2 = -x12*dnorm(2,j) + x22*dnorm(1,j);

ds = abs(d);
dtheta = atan2(ds*al, ds^2+t1*t2);

aa = (-dtheta*ds + al + t1*log(r1) - t2*log(r2)) / (pi*2);

if d < 0
    dtheta = -dtheta;
end
if i == j
    bb = 0.5;
else
    bb = -dtheta / (pi*2);
end


return
end


