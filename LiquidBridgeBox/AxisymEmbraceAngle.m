function [alpha] = AxisymEmbraceAngle(R,r)
% Calculate embracing angle

if R > 0
	alpha = asin(r/R);
else
	alpha = 0;
end

return
end
