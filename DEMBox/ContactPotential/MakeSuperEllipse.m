function [shape] = MakeSuperEllipse(a,b,p,q, xc,yc,rotang)

shape = struct();

shape.type = 2;

shape.a = a;
shape.b = b;
shape.p = p;
shape.q = q;

shape.xc = xc;
shape.yc = yc;
shape.rotang = rotang;

return
end

