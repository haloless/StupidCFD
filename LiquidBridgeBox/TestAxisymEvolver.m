
clear all;

sigma = 1.0;
R1 = 1.0;
R2 = 1.0;
H = 0.0;
theta1 = deg2rad(0);
theta2 = deg2rad(0);
% V = SphereVolume(R1)*0.2;
V = 1.0e-7;

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










