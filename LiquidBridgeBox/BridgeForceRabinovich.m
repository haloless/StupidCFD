function [F] = BridgeForceRabinovich(R,H,theta,V,sigma)

if H > 0
	h2d = 1 / (-1 + sqrt(1+2*V/(pi*R*H^2)));
else
	h2d = 0;
end

F = 2*pi*R*sigma*cos(theta) / (1+h2d);

return
end

