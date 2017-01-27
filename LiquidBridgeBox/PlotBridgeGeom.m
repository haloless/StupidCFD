function [] = PlotBridgeGeom(bridge)

R1 = bridge.R1;
R2 = bridge.R2;
H = bridge.H;
X1 = bridge.X1;
X2 = bridge.X2;

% draw sphere 1 
PlotCircle([X1,0],R1, 'EdgeColor','k');

if R2 > 0
	% draw sphere 2
	PlotCircle([X2,0],R2, 'EdgeColor','k');
else
	% draw straight wall
	line([X2,X2], [0,R1*1.5], 'Color','k');
end

axis equal;
axis([X1 X2 0 R1*1.25]);


return
end
