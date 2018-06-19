function [R] = DerjaguinRadius(R1,R2)
% Average different radii
% For wall, R2->infty

if R2 > 0
	R = 2*R1*R2/(R1+R2);
else
	R = 2*R1;
end

return
end


