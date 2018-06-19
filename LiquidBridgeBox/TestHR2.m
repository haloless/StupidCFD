%%
%% 
%%

clear all;

sigma = 1.0;
R1 = 1;
% R1 = 2;
R2 = 1;
% R2 = 2;
% R2 = 3;
% R2 = -1;
theta1 = deg2rad(0);
% theta1 = deg2rad(15);
% theta1 = deg2rad(30);
% theta1 = deg2rad(60);
% theta1 = deg2rad(90);
% theta1 = deg2rad(120);
% theta1 = deg2rad(150);
theta2 = deg2rad(0);
% theta2 = deg2rad(15);
% theta2 = deg2rad(30);
% theta2 = deg2rad(60);
% theta2 = deg2rad(90);
% theta2 = deg2rad(100);
% theta2 = deg2rad(120);
% theta2 = deg2rad(150);

V1 = 4/3*pi*R1^3;
% V = V1 * 0.0001;
% V = V1 * 0.0005;
% V = V1 * 0.001;
% V = V1 * 0.005;
% V = V1 * 0.01;
% V = V1 * 0.05;
% V = V1 * 0.1;
% V = V1 * 0.2;
% V = V1 * 0.3;
% V = V1 * 0.4;

% V = 0.01;
% V = 0.04;
% V = 0.05;
% V = 0.1;
V = 0.25;

H = 0.0;
% H = 0.05;
% H = 0.1;
% H = 0.2;
% H = 0.25;

% Hrup1 = BridgeRuptureMikami(R1,R2,theta1,V)
% Hrup2 = BridgeRuptureLian(theta1,V)
% F = BridgeForceHR2(R1,R2,H,theta1,theta2,V,sigma);

% F1 = BridgeForceWillet(R1,H,theta1,V,sigma);

bridge = MakeBridge(R1,R2,H,theta1,theta2,V,sigma);
% Fref = AxisymEvolverDirectForce(bridge,np,[],[]);

%
X1 = bridge.X1;
X2 = bridge.X2;

% draw spheres
hfig = figure;
hold on;
PlotBridgeGeom(bridge);
hold off;

%
% np = 31;
np = 51;
% np = 81;
% np = 101;
% np = 201;
[rp,xp] = AxisymEvolverGuessInit(bridge,np);

if 0
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
Fopt = AxisymEvolverEvalForce(bridge, np,rp,xp,pres);


if 1
	% draw optimized shape
	figure(hfig);
	hold on;
	plot(xp,rp,'.-r');
	hold off;
	drawnow;
end 


[F,res] = BridgeForceHR2(R1,R2,H,theta1,theta2,V,sigma);
if 1
    figure(hfig);
    hold on;
    if res.rout > 0
        rectangle('Position',[res.xcirc-res.rout,res.ycirc-res.rout,res.rout*2,res.rout*2], 'Curvature',[1,1], ...
        'EdgeColor','b');
    else
        rectangle('Position',[res.xcirc+res.rout,res.ycirc+res.rout,-res.rout*2,-res.rout*2], 'Curvature',[1,1], ...
        'EdgeColor','b');
    end
    plot(res.xcirc,res.ycirc,'xb');
    hold off;
end





% F
% F1
% Fref


