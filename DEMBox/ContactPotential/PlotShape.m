
function [h] = PlotShape(shape,plotcell)

% nseg = 360;

% theta = (0:nseg) .* (2*pi/nseg);

% ea = shape.a;
% eb = shape.b;
% xc = shape.xc;
% yc = shape.yc;
% rotang = shape.rotang;

% xs = ea .* cos(theta);
% ys = eb .* sin(theta);
% pos = [xs;ys];


% rotmat = RotationMatrix(rotang);

% pos = rotmat * pos;
% xs = pos(1,:) + xc;
% ys = pos(2,:) + yc;

% plot(xs,ys,plotcell{:});

% range = [shape.xc-shape.a,shape.xc+shape.a, shape.yc-shape.b, shape.yc+shape.b];

if shape.type >= 0
h = ezplot(@(x,y) ShapePotential(shape,x,y));
else
h = ezplot(@(x,y) WallPotential(shape,x,y));
end
set(h, 'Color',plotcell{1});



return
end

