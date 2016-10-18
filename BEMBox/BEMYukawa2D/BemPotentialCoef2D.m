

function [ aa,bb ] = BemPotentialCoef2D(xin,yin,jelem)

BemMeshGlobals;

% element
al = dlen(jelem);
j1 = node(1,jelem);
j2 = node(2,jelem);

% source point
x11 = y(1,j1) - xin;
x21 = y(2,j1) - yin;
x12 = y(1,j2) - xin;
x22 = y(2,j2) - yin;

r1 = sqrt(x11^2 + x21^2);
r2 = sqrt(x12^2 + x22^2);
d = x11*dnorm(1,jelem) + x21*dnorm(2,jelem);
t1 = x11*dtang(1,jelem) + x21*dtang(2,jelem);
t2 = x12*dtang(1,jelem) + x22*dtang(2,jelem);

ds = abs(d);
theta1 = atan2(t1,ds);
theta2 = atan2(t2,ds);
dtheta = theta2 - theta1;

aa = (-dtheta*ds + al + t1*log(r1) - t2*log(r2)) / (pi*2);
if d < 0
    dtheta = -dtheta;
end
bb = -dtheta / (pi*2);

return
end




