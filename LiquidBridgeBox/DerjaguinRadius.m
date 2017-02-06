function [R] = DerjaguinRadius(R1,R2)

if R2 > 0
	R = 2*R1*R2/(R1+R2);
else
	R = 2*R1;
end

return
end


