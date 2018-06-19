function [] = PlotBridgeGeom(bridge)
% Plot the fixed spheres (walls).
% Remember to activate the figure and hold on.

% parameters
R1 = bridge.R1;
R2 = bridge.R2;
H = bridge.H;
X1 = bridge.X1;
X2 = bridge.X2;

% draw sphere 1 
% PlotCircle([X1,0],R1, 'EdgeColor','k');
% ts = linspace(0,pi/2,361);
ts = linspace(0,pi/4,361);
plot(R1*cos(ts),R1*sin(ts),'k-');

if R2 > 0
	% draw sphere 2
	% PlotCircle([X2,0],R2, 'EdgeColor','k');
    % ts = linspace(pi/2,pi,361);
    ts = linspace(pi/4*3,pi,361);
    plot(X2+R2*cos(ts),R2*sin(ts),'k-');
else
	% draw straight wall
	line([X2,X2], [0,R1*1.5], 'Color','k');
end

axis equal;
axis([X1 X2 0 R1*1.25]);


return
end
