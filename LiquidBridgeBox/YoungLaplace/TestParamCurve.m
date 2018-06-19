
clear;

sigma = 1.0;
R1 = 1.0;
R2 = 1.0;
R2 = -1.0;

theta1 = deg2rad(0);
% theta1 = deg2rad(20);
% theta1 = deg2rad(30);
% theta1 = deg2rad(40);
% theta1 = deg2rad(60);
% theta2 = deg2rad(0);
theta2 = deg2rad(10);
% theta2 = deg2rad(20);
% theta2 = deg2rad(30);
% theta2 = deg2rad(40);
% theta2 = deg2rad(60);

% V = 0.05;
% V = 0.0348;
% V = 0.008;
% V = 0.3^3;
% V = 0.0001;
% V = 0.001;
% V = 0.01;
% V = 0.04;
V = 0.06;
% V = 0.1;

% H = 0.2;
H = 0.3;
% H = 0.4;
H = 0.376;
% H = 0.273;
% H = 0.265;
% H = 0.3072;
% H = 0.42;
% H = 0.0365;
% H = 0.0765;
% H = 0.170;
% H = 0.271;
% H = 0.370;
% H = 0.67;
% H = 0.505;
% H = 0.4302;

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
if 1
    RunEvolver;
    alpha1 = asin(rp(1)/R1);
    alpha2 = asin(rp(end)/R2);
end


if 0
    % alpha1 = asin(rp(1)/R1);
    % alpha2 = asin(rp(end)/R2);
    % pres = pres;
    
    [~,guess] = BridgeForceHR2(R1,R2,H,theta1,theta2,V,sigma);
    alpha1 = guess.alpha1;
    alpha2 = guess.alpha2;
    pres = guess.pres
    
    tic;
    [alpha1,alpha2,pres] = ParamSolvePres(bridge, alpha1,alpha2,pres);
    toc;
    
    pres
end

if 0
    M = -pres / 2;
    C = R1*sin(alpha1)*sin(alpha1+theta1) + (R1^2)*M*sin(alpha1)^2;
    
    phi1 = -(pi/2 - alpha1 - theta1);
    phi2 = (pi/2 - alpha2 - theta2);
    
    phis = linspace(phi1,phi2, 51);
    phis = linspace(phi1,phi1*1.0015,51);
    
    [xs,ys] = ParamCurve(phis, phi1, C,M);
    % remember to shift x
    xs = xs + R1*cos(alpha1);
    
    figure(hfig);
    hold on;
    plot(xs,ys,'.-');
    hold off;
end

if 0
    vol = ParamVolume(phi1,phi2,C,M);
    vol1 = SphereCapVolume(R1*sin(alpha1),R1*(1-cos(alpha1)));
    vol2 = SphereCapVolume(R2*sin(alpha2),R2*(1-cos(alpha2)));
    vol-vol1-vol2
end


if 0
    [ts,zs] = YLIntegrateBridge(bridge, alpha1,pres);
    
    figure(hfig);
    hold on;
    plot(ts,zs(:,1),'k+-');
    hold off;

end

if 0
    tic;
    sol = YLSolveBVP(bridge);
    toc;
    
    sol.ok
    sol.parameters
    
    figure(hfig);
    hold on;
    tmps = linspace(0,1,81);
    if R2 > 0
        L = R1*(1-cos(sol.parameters(1))) + H + R2*(1-cos(sol.parameters(2)));
    else
        L = R1*(1-cos(sol.parameters(1))) + H;
    end
    xs = R1*cos(sol.parameters(1)) + L.*tmps;
    zs = deval(sol, tmps);
    plot(xs, zs(1,:), '.-b');
    hold off;
end

if 1
    tic;
    if 1
        sol = YLShootBVP(bridge);
    else
        guess = [alpha1,alpha2,pres];
        if R2 <= 0
            guess(2) = rp(end);
        end
        sol = YLShootBVP(bridge, guess);
    end
    toc;
    
    sol
    
    [ts,zs] = YLIntegrateBridge(bridge, sol.alpha1,sol.pres);
    
    figure(hfig);
    hold on;
    plot(ts,zs(:,1),'k.-');
    hold off;
end



