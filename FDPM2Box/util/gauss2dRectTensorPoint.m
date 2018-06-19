function [gn,gp,gw] = gauss2dRectTensorPoint(gnx, gny)
%gauss2dRectTensorPoint: 2D tensor product of 1D gauss  


if nargin == 1
    gny = gnx;
end

% create 1d points
[gpx,gwx] = gauss1dPoint(gnx);
[gpy,gwy] = gauss1dPoint(gny);

% tensor product
gx = zeros(gnx,gny);
gy = zeros(gnx,gny);
gw = zeros(gnx,gny);

for j = 1:gny
for i = 1:gnx
    gx(i,j) = gpx(i);
    gy(i,j) = gpy(j);
    gw(i,j) = gwx(i) * gwy(j);
end
end

gn = gnx * gny;
gp = [gx(:),gy(:)];
gw = gw(:);


return
end
