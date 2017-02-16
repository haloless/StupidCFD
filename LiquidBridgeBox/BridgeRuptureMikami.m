function [hrup] = BridgeRuptureMikami(R1,R2,theta,V)

if R2 > 0
	% sphere-sphere	
	hrup = (0.99+0.62*theta) * V^(0.34);
else
	% sphere-wall
	hrup = (0.95+0.22*theta) * V^(0.32);
end

return
end


