function [shape] = MakeEllipse(a,b, xc,yc,rotang)

shape = struct();

shape.type = 1;

shape.a = a;
shape.b = b;

shape.xc = xc;
shape.yc = yc;
shape.rotang = rotang;

return
end

