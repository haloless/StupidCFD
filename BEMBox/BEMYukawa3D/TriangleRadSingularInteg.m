function [qint] = TriangleRadSingularInteg(qfunc,x0,x1,x2, nga,wga,xga, ngr,wgr,xgr)

h = 1/sqrt(2);

qint = 0;
for i = 1:nga
	% angle
	phi = pi/4 * xga(i);
	cphi = cos(phi);
	sphi = sin(phi);
	cang = cos(phi+pi/4);
	sang = sin(phi+pi/4);
	wphi = wga(i);
	
	for j = 1:ngr
		% radius
		xi = (xgr(j)+1)/2;
		rad = xi * h/cphi;
		wrad = wgr(j);
		
		e1 = rad * cang;
		e2 = rad * sang;
		qpos = x0 + e1*(x1-x0) + e2*(x2-x0);
		
		qval = qfunc(qpos(1),qpos(2),qpos(3));
		
		qint = qint + qval*rad * wphi*wrad * (pi/4) * (0.5*h/cphi);
	end
end

% patch area
% note coef 0.5 is already treated in preceding integration
aint = cross(x1-x0,x2-x0);
% aint = 0.5 * aint(3);
% aint = aint(3);
aint = norm(aint);

qint = qint * aint;


return
end

