function [phi] = ShapePotential(shape,x,y)

% common operation
% shift to center
xp = x - shape.xc;
yp = y - shape.yc;
% back-transform to body coordinate
rmat = RotationMatrix(-shape.rotang);
pp = rmat * [xp';yp'];
xp = pp(1,:)';
yp = pp(2,:)';

switch shape.type
	case 1
		phi = xp.^2/shape.a^2 + yp.^2/shape.b^2 - 1;
	case 2
		phi = (xp./shape.a).^shape.p + (yp./shape.b).^shape.q - 1;
	otherwise
		error(['Unknown shape.type=',int2str(shape.type)]);
end


return
end

