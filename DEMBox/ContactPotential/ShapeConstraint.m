function [c,ceq] = ShapeConstraint(shape,pos)
c = [];
ceq = EllipsePotential(shape,pos(1),pos(2));
return
end


