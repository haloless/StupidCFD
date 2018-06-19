
clear;

sigma = 1.0;
R1 = 1.0;
% R2 = 1.0;
R2 = -1.0;
H = 1.878935e-1;
% H = 0;
% theta1 = deg2rad(80.6);
% theta2 = deg2rad(80.6);
theta1 = deg2rad(60);
theta2 = deg2rad(100);
V = 0.08;
% V = 0.01;
% V = 0.8;
% V = 4/3*pi*R1^3 * 0.4;
% V = 0.15;

thetam = acos(0.5*(cos(theta1)+cos(theta2)));

bridge = MakeBridge(R1,R2,H,theta1,theta2,V,sigma);

%
X1 = bridge.X1;
X2 = bridge.X2;

% draw spheres
hfig = figure;
hold on;
PlotBridgeGeom(bridge);
hold off;

%
np = 51;
[rp,xp] = AxisymEvolverGuessInit(bridge,np);

if 1
	% draw initial shape
	figure(hfig);
	hold on;
	plot(xp,rp,'-b');
	hold off;
	drawnow;
end
% return

% optimize
[rp,xp,pres] = AxisymEvolver(bridge,np,rp,xp);

if 1
	% draw optimized shape
	figure(hfig);
	hold on;
	plot(xp,rp,'.-r');
	hold off;
	drawnow;
end 

% force
[F1,~,~,F2,~,~] = AxisymEvolverEvalForce(bridge, np,rp,xp,pres)










