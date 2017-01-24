function [qint] = TriangleGaussInteg(qfunc,x0,x1,x2, ng,wg,xg,yg)

qint = 0;

qx = (1-xg-yg).*x0(1) + xg.*x1(1) + yg.*x2(1);
qy = (1-xg-yg).*x0(2) + xg.*x1(2) + yg.*x2(2);
qz = (1-xg-yg).*x0(3) + xg.*x1(3) + yg.*x2(3);
qval = qfunc(qx,qy,qz);
qint = sum(qval.*wg);

% for ig = 1:ng
	% qpos = x0 + xg(ig)*(x1-x0) + yg(ig)*(x2-x0);
	% qval = qfunc(qpos(1),qpos(2),qpos(3));
	% qint = qint + qval*wg(ig);
% end

% triangle area
aint = cross(x1-x0,x2-x0);
aint = norm(aint);
aint = 0.5 * aint;

qint = qint * aint;


return
end


