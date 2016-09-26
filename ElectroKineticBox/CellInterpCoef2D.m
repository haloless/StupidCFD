
function [ is,js,ws ] = CellInterpCoef2D(xint,yint, nx,ny,dx,dy,xlo,ylo)

xbase = xlo + dx/2;
ybase = ylo + dy/2;

i0 = floor((xint-xbase)/dx) + 1;
j0 = floor((yint-ybase)/dy) + 1;

x0 = (i0-1)*dx + xbase;
y0 = (j0-1)*dy + ybase;

rx = (xint-x0) / dx;
ry = (yint-y0) / dy;

is = [ i0, i0+1, i0,   i0+1 ];
js = [ j0, j0,   j0+1, j0+1 ];
ws = [(1-rx)*(1-ry), rx*(1-ry), (1-rx)*ry, rx*ry ];


return
end

