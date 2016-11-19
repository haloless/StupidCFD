
function [ ws,is,js ] = GridInterpCoef2D(xint,yint, xlo,ylo,dx,dy)

xbase = xlo + dx/2;
ybase = ylo + dy/2;

i0 = floor((xint-xbase) / dx) + 1;
j0 = floor((yint-ybase) / dy) + 1;

x0 = (i0-1)*dx + xbase;
y0 = (j0-1)*dy + ybase;

rx = (xint-x0) / dx;
ry = (yint-y0) / dy;

is(1) = i0;
js(1) = j0;
ws(1) = (1-rx) * (1-ry);
is(2) = i0 + 1;
js(2) = j0;
ws(2) = rx * (1-ry);
is(3) = i0;
js(3) = j0 + 1;
ws(3) = (1-rx) * ry;
is(4) = i0 + 1;
js(4) = j0 + 1;
ws(4) = rx * ry;



return
end
