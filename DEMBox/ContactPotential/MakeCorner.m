function [corner] = MakeCorner(xw,xa,xb)

corner = struct();

corner.type = -2;

corner.xw = xw;
corner.xa = xa;
corner.xb = xb;

t1 = normalize(xw - xa);
n1 = [-t1(2), t1(1)]';
corner.t1 = t1;
corner.n1 = n1;
corner.d1 = dot(n1,xw);

t2 = normalize(xb - xw);
n2 = [-t2(2), t2(1)]';
corner.t2 = t2;
corner.n2 = n2;
corner.d2 = dot(n2,xw);

nw = n1 + n2;
nw = nw ./ norm(nw);
corner.nw = nw;

return
end

function [vn] = normalize(v)
    vn = v ./ norm(v);
return
end

