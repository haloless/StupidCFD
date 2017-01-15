function [gw,hw] = CornerDeriv(corner,xc,yc,dh)

dc = CornerPotential(corner,xc,   yc   );
dw = CornerPotential(corner,xc-dh,yc   );
de = CornerPotential(corner,xc+dh,yc   );
ds = CornerPotential(corner,xc   ,yc-dh);
dn = CornerPotential(corner,xc   ,yc+dh);

wx = (de-dw) / (dh*2);
wy = (dn-ds) / (dh*2);
wxx = (de+dw-dc*2) / (dh*dh);
wyy = (dn+ds-dc*2) / (dh*dh);

d00 = CornerPotential(corner,xc-dh,yc-dh);
d01 = CornerPotential(corner,xc+dh,yc-dh);
d10 = CornerPotential(corner,xc-dh,yc+dh);
d11 = CornerPotential(corner,xc+dh,yc+dh);

wxy = ((d11-d10)/(dh*2) - (d01-d00)/(dh*2)) / (dh*2);

gw = [ wx; wy ];
hw = [ wxx, wxy; wxy, wyy ];

return
end


