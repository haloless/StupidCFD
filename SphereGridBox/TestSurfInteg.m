function [sint] = TestSurfInteg(nref)

% clear all;

nodes = IcosahedralPoints(nref);
tri = SpherePointsToGrid(nodes(1,:),nodes(2,:),nodes(3,:));
ntri = size(tri,1);

% pos1 = nodes(:,tri(1,1));
% pos2 = nodes(:,tri(1,2));
% pos3 = nodes(:,tri(1,3));
% pos4 = (pos1+pos2)/2; pos4 = pos4 ./ norm(pos4);
% pos5 = (pos2+pos3)/2; pos5 = pos5 ./ norm(pos5);
% pos6 = (pos3+pos1)/2; pos6 = pos6 ./ norm(pos6);
% pos = [pos1,pos2,pos3,pos4,pos5,pos6];

area = 0;
for iface = 1:ntri
	v1 = nodes(:,tri(iface,1));
	v2 = nodes(:,tri(iface,2));
	v3 = nodes(:,tri(iface,3));
	area = area + 0.5*norm(cross(v2-v1,v3-v1));
end
area

if (1)
	figure;
	trimesh(tri,nodes(1,:),nodes(2,:),nodes(3,:));
	axis equal;
	hold on;
	% plot3(pos1(1),pos1(2),pos1(3),'x');
	% plot3(pos2(1),pos2(2),pos2(3),'x');
	% plot3(pos3(1),pos3(2),pos3(3),'x');
	% plot3(pos(1,:),pos(2,:),pos(3,:),'x');
	% for i = 1:6
		% text(pos(1,i),pos(2,i),pos(3,i), [int2str(i)]);
	% end
	hold off;
end

ng = 7;
[cg,ug,vg] = TriangleGaussRule(ng);

aint = 0;

for iface = 1:ntri
	pos1 = nodes(:,tri(iface,1));
	pos2 = nodes(:,tri(iface,2));
	pos3 = nodes(:,tri(iface,3));
	% map midpoint to sphere surface
	pos4 = (pos1+pos2)/2; pos4 = pos4 ./ norm(pos4);
	pos5 = (pos2+pos3)/2; pos5 = pos5 ./ norm(pos5);
	pos6 = (pos3+pos1)/2; pos6 = pos6 ./ norm(pos6);
	% six point element
	pos = [pos1,pos2,pos3,pos4,pos5,pos6];
	
	sint = 0;
	for ig = 1:ng
		c = cg(ig);
		u = ug(ig);
		v = vg(ig);
		
		[N,Nu,Nv] = ShapeFunc(u,v);
		
		posg = pos * N;
		eu = pos * Nu;
		ev = pos * Nv;
		
		ds = 0.5 * norm(cross(eu,ev));
		sval = 1.0;
		sint = sint + sval*ds*c;
	end
	aint = aint + sint;
end
aint

return
end

function [N,Nu,Nv] = ShapeFunc(u,v)
	w = 1-u-v;
	N = [ w*(2*w-1); u*(2*u-1); v*(2*v-1); 4*w*u; 4*u*v; 4*v*w ];
	Nu = [ -4*w+1; 4*u-1; 0; 4*w-4*u; 4*v; -4*v ];
	Nv = [ -4*w+1; 0; 4*v-1; -4*u; 4*u; 4*w-4*v ];
	return
end




function [wg,xg,yg] = TriangleGaussRule(ng)

if ng == 7
	wg = [...
	0.225000000000000;
	0.125939180544827;
	0.125939180544827;
	0.125939180544827;
	0.132394152788506;
	0.132394152788506;
	0.132394152788506];
	
	xg = [...
	0.333333333333333;
	0.797426985353087;
	0.101286507323456;
	0.101286507323456;
	0.470142064105115;
	0.470142064105115;
	0.059715871789770];
	
	yg = [...
	0.333333333333333;
	0.101286507323456;
	0.797426985353087;
	0.101286507323456;
	0.470142064105115;
	0.059715871789770;
	0.470142064105115];
	
else
	error(['Unsupported Gaussian Quadrature NG=',int2str(ng)]);
end


return
end
