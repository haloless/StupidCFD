
clear all;


% shape1 = MakeSuperEllipse(3.0,1.5, 2,2, -0.0,0.0, 40/180*pi);

% shape1 = MakeSuperEllipse(3.0,1.5, 4,4, -0.0,0.0, 40/180*pi);
% shape1 = MakeSuperEllipse(3.0,1.5, 6,6, -0.0,0.0, 40/180*pi);
% shape1 = MakeSuperEllipse(3.0,1.5, 2,2, -0.0,0.0, 0/180*pi);
shape1 = MakeSuperEllipse(3.0,1.5, 6,6, -0.0,0.0, 45/180*pi);
% shape1 = MakeSuperEllipse(3.0,1.5, 6,6, -0.0,0.0, 0/180*pi);

% wall2 = MakeWall(1.0,-1.0, -1.0,1.0);
% wall2 = MakeWall(1.0,-1.0, 0.0,1.0);
% wall2 = MakeWall(1.0,-1.0, -1.0,0.0);
% wall2 = MakeWall(1.0,-2.0, -1.0,2.0);
wall2 = MakeWall(1.0,-1.0, -1.0,1.0);

hfig = figure;
hold on;
PlotShape(shape1, {'b'});
PlotShape(wall2, {'r'});
plot(wall2.xwall(1),wall2.xwall(2),'.r');
quiver(wall2.xwall(1),wall2.xwall(2),wall2.nwall(1),wall2.nwall(2),'r');
axis([-6 6 -3 3]);
axis equal;
hold off;


if 0
	figure(hfig);
	hold on;
	for x = linspace(-2,6,50)
	for y = linspace(-3,3,40)
		phi1 = ShapePotential(shape1,x,y);
		phi2 = WallPotential(wall2,x,y);
		if phi1 <= 0
			plot(x,y,'xb');
		end
		if phi2 <= 0
			plot(x,y,'or');
		end
	end
	end
	hold off;
end

fun1 = @(xy) ShapePotential(shape1,xy(1),xy(2));
fun2 = @(xy) WallPotential(wall2,xy(1),xy(2));

con1 = @(xy) deal([], ShapePotential(shape1,xy(1),xy(2)));
con2 = @(xy) deal([], WallPotential(shape2,xy(1),xy(2)));

% search from mid point
xmid = 0.5 * (shape1.xc + wall2.xwall(1));
ymid = 0.5 * (shape1.yc + wall2.xwall(2));
guess = [xmid;ymid];

options = optimoptions('fmincon', ...
'Display','iter', 'Algorithm','sqp', 'MaxFunEvals',2000);

[sol1,phi1] = fmincon(fun2,guess,[],[],[],[],[],[], con1, options);
% [sol2,phi2] = fmincon(fun1,guess,[],[],[],[],[],[], con2, options);

% sol3 = SolveContactPotential(shape1,shape2);
% sol4 = SolveContactPotential(shape2,shape1);

sol2 = SolveContactWall(shape1,wall2);



if 1
	figure(hfig);
	hold on;
	
	plot(sol1(1),sol1(2), 'xb');
	plot(sol2(1),sol2(2), 'ob');
	
	hold off;
end





