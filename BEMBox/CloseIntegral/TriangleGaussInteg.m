function [qint] = TriangleGaussInteg(qfunc,x0,x1,x2, ng,wg,xg,yg)

qint = 0;

for ig = 1:ng
	qpos = x0 + xg(ig)*(x1-x0) + yg(ig)*(x2-x0);
	qval = qfunc(qpos(1),qpos(2));
	qint = qint + qval*wg(ig);
end

% triangle area
aint = cross([x1-x0;0],[x2-x0;0]);
aint = 0.5 * aint(3);

qint = qint * aint;


return
end


