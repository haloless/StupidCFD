
clear all;

% shape1 = MakeEllipse(3.0,1.5, -3.0,0.0, -15/180*pi);
% shape2 = MakeEllipse(4.0,2.5, 3.0,1,  30/180*pi);
% shape2 = MakeEllipse(4.0,2.5, 3.5,1,  30/180*pi);

% shape1 = MakeSuperEllipse(3.0,1.5, 2,2, -3.0,0.0, -15/180*pi);
% shape2 = MakeSuperEllipse(4.0,2.5, 6,6, 3.5,1,  40/180*pi);
% shape2 = MakeSuperEllipse(4.0,2.5, 6,6, 3.8,1,  30/180*pi);
% shape2 = MakeSuperEllipse(4.0,2.5, 6,6, 4.8,1,  30/180*pi);

% shape1 = MakeSuperEllipse(3.0,1.5, 2,2, -3.0,0.0, 90/180*pi);
shape1 = MakeSuperEllipse(3.0,1.5, 6,6, -3.0,0.0, 90/180*pi);
% shape2 = MakeSuperEllipse(4.0,2.5, 6,6, 2.5,1,  30/180*pi);
% shape2 = MakeSuperEllipse(4.0,2.5, 4,4, 3.0,1,  30/180*pi);
shape2 = MakeSuperEllipse(4.0,2.5, 4,4, 3.0,0,  30/180*pi);
% shape2 = MakeSuperEllipse(4.0,2.5, 2,2, 2.0,1,  0/180*pi);

% shape1 = MakeSuperEllipse(3.0,1.5, 6,6, -3.0,0.0, 90/180*pi);
% shape2 = MakeSuperEllipse(4.0,2.5, 6,6, 2.5,1,  0/180*pi);
% shape2 = MakeSuperEllipse(4.0,2.5, 2,2, 2.0,1,  0/180*pi);

% shape1 = MakeSuperEllipse(3.0,1.5, 6,6, -3.0,0.0, 0/180*pi);
% shape2 = MakeSuperEllipse(4.0,2.5, 6,6, 2.5,0.0,  0/180*pi);


hfig = figure;
hold on;
hs1 = PlotShape(shape1, {'b'});
hs2 = PlotShape(shape2, {'r'});
axis equal;
axis([-6 6 -5 5]);
hold off;

[node1,elem1] = MakeShapeGrid(shape1);
[node2,elem2] = MakeShapeGrid(shape2);

figure(hfig);
hold on;
plot(node1(1,:),node1(2,:),'.b');
plot(node2(1,:),node2(2,:),'.r');
hold off;

if 0
	% build
    sdf1 = MakeGridSDF(shape1,node1,elem1);
    sdf2 = MakeGridSDF(shape2,node2,elem2);
	save('cache/shape1.mat','sdf1');
	save('cache/shape2.mat','sdf2');
else
	% load
	load('cache/shape1.mat','sdf1');
	load('cache/shape2.mat','sdf2');
end

shape2.xc = shape2.xc - 0.5;
sdf2.xg = sdf2.xg - 0.5;
sdf2.xmin = sdf2.xmin - 0.5;
sdf2.xmax = sdf2.xmax - 0.5;

if 1
    % do not converge without reinit
    sdf1.phig = ImplicitFuncReinit(sdf1.phig,sdf1.nx,sdf1.ny,sdf1.dh,sdf1.dh, 2.0);
    sdf2.phig = ImplicitFuncReinit(sdf2.phig,sdf2.nx,sdf2.ny,sdf2.dh,sdf2.dh, 2.0);
end

figure(hfig);
hold on;
PlotSDF(sdf1, [0,0],'g','ShowText','on'); 
PlotSDF(sdf2, [0,0],'k','ShowText','on'); 
% PlotSDF(sdf1, 'g','ShowText','on'); 
% PlotSDF(sdf2, 'k','ShowText','on'); 
hold off;

% return

xcont = SolveContactSDF(sdf1,sdf2, [(shape1.xc+shape2.xc)/2;(shape1.yc+shape2.yc)/2]);

figure(hfig);
hold on;
plot(xcont(1),xcont(2),'rx');
hold off;


return



if 0
	figure(hfig);
	hold on;
	for x = linspace(-2,6,50)
	for y = linspace(-3,3,40)
		phi1 = ShapePotential(shape1,x,y);
		phi2 = ShapePotential(shape2,x,y);
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
fun2 = @(xy) ShapePotential(shape2,xy(1),xy(2));

con1 = @(xy) deal([], ShapePotential(shape1,xy(1),xy(2)));
con2 = @(xy) deal([], ShapePotential(shape2,xy(1),xy(2)));

% search from mid point
xmid = 0.5 * (shape1.xc + shape2.xc);
ymid = 0.5 * (shape1.yc + shape2.yc);

guess = [xmid;ymid];
options = optimoptions('fmincon', 'Display','iter', 'Algorithm','sqp');

[sol1,phi1] = fmincon(fun2,guess,[],[],[],[],[],[], con1, options);
[sol2,phi2] = fmincon(fun1,guess,[],[],[],[],[],[], con2, options);

sol3 = SolveContactPotential(shape1,shape2);
sol4 = SolveContactPotential(shape2,shape1);

if 1
	figure(hfig);
	hold on;
	
	plot(sol1(1),sol1(2), 'xb');
	plot(sol2(1),sol2(2), 'xr');
	
	plot(sol3(1),sol3(2), 'ob');
	plot(sol4(1),sol4(2), 'or');
	
	% text(sol1(1),sol1(2), num2str(phi1));
	% text(sol2(1),sol2(2), num2str(phi2));
	
	hold off;
end





