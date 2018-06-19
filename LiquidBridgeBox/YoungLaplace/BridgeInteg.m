function [res] = BridgeInteg(bridge, alpha1, pres)
	
	R1 = bridge.R1;
	R2 = bridge.R2;
	H = bridge.H;
	theta1 = bridge.theta1;
	theta2 = bridge.theta2;
	X1 = bridge.X1;
	X2 = bridge.X2;
	
	if R2 <= 0
		R2 = 0;
	end
	
	%
	M = -pres / 2;
	C = R1*sin(alpha1)*sin(alpha1+theta1) + M*R1^2*sin(alpha1)^2;
	
	% ODE function
	function [dz] = yl_ode(t,z)
		y = z(1);
		dy = z(2);
		q = sqrt(1+dy.^2);
		ddy = q.^2 ./ y + 2*M * q.^3;
		
		dz = zeros(2,1);
		dz(1) = dy;
		dz(2) = ddy;
	return
	end
	
	% Event function
	function [val,isterm,dir] = yl_event(t,z)
		x = t;
		y = z(1);
		dy = z(2);
		
		if R2 > 0
			r = sqrt((x-X2).^2 + y.^2);
			dr = ((x-X2) + y.*dy) ./ r;
			val = (r-R2) .* dr;
		else
			val = x - X2;
		end
		
		if (abs(dy) > 1.0e3)
			val = 0;
		end
		
		isterm = 1;
		dir = 0;
		
		
	return
	end
	
	
	%
	phi1 = -(pi/2 - alpha1 - theta1);
	x1 = R1 * cos(alpha1);
	y1 = R1 * sin(alpha1);
	dy1 = tan(phi1);
	
	%
	% options = odeset('RelTol',1.0e-6, 'AbsTol',1.0e-6);
	options = odeset('RelTol',1.0e-10, 'AbsTol',1.0e-10, 'Events',@yl_event);
	[xs,zs] = ode45(@yl_ode, [x1;x1+R1+H+R2], [y1;dy1], options);
	
	% end point
	x2 = xs(end);
	y2 = zs(end,1);
	dy2 = zs(end,2);
	
	%
	alpha2 = atan2(y2, X2-x2);
	phi2 = atan(dy2);
	ca2 = pi/2 - alpha2 - phi2;
	
	% if end point on sphere 2
	touch = 0;
	dist2 = sqrt((x2-X2)^2+y2^2);
	if abs(dist2-R2) < 1.0e-6
		touch = 1;
		if abs(ca2-theta2) < 1.0e-4
			touch = 2;
		end
	end
	
	%
	res.xs = xs;
	res.ys = zs(:,1);
	res.dys = zs(:,2);
	%
	res.x2 = x2;
	res.y2 = y2;
	res.dy2 = dy2;
	res.alpha2 = alpha2;
	res.phi2 = phi2;
	res.ca2 = ca2;
	%
	res.touch = touch;
	
return
end


