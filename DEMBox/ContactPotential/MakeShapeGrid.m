function [node,elem] = MakeShapeGrid(shape)

if shape.type~=2
    error('type');
end

n = 36;
theta = 2*pi/n * (0:n-1);
xx = cos(theta);
yy = sin(theta);
sx = sign(xx);
sy = sign(yy);
px = shape.a * (xx.^2).^(1/shape.p) .* sx;
py = shape.b * (yy.^2).^(1/shape.q) .* sy;

node = [px; py];

rmat = RotationMatrix(shape.rotang);
node = rmat * node;
node(1,:) = node(1,:) + shape.xc;
node(2,:) = node(2,:) + shape.yc;

elem = [1:n; 2:n+1]';
elem(n,2) = 1;


return
end

