function [cc] = triCircumCenter(ax,ay, bx,by, cx,cy)
%triCircumCenter: Compute circumcenter of triangle
% CC is a row vector

bcy = by - cy;
cay = cy - ay;
aby = ay - by;

cbx = cx - bx;
acx = ax - cx;
bax = bx - ax;

a2 = ax.^2 + ay.^2;
b2 = bx.^2 + by.^2;
c2 = cx.^2 + cy.^2;

dd = 2 .* (ax.*bcy + bx.*cay + cx.*aby);

uu = (a2.*bcy + b2.*cay + c2.*aby) ./ dd;
vv = (a2.*cbx + b2.*acx + c2.*bax) ./ dd;

cc = [uu,vv];


return
end

