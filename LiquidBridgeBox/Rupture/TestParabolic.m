
clear;

sigma = 1.0;

R1 = 1.0;
R2 = 1.0;
% R2 = -1.0;

Rm = DerjaguinRadius(R1,R2);

% theta1 = deg2rad(0);
% theta1 = deg2rad(10);
% theta1 = deg2rad(15);
% theta1 = deg2rad(20);
theta1 = deg2rad(30);
% theta1 = deg2rad(40);
% theta1 = deg2rad(60);
% theta2 = deg2rad(0);
% theta2 = deg2rad(10);
% theta2 = deg2rad(15);
% theta2 = deg2rad(20);
theta2 = deg2rad(30);
% theta2 = deg2rad(40);
% theta2 = deg2rad(50);
% theta2 = deg2rad(60);

% V = 0.001;
% V = 0.003;
% V = 0.006;
V = 0.01;
% V = 0.03;
% V = 0.06;

% Hs = 0.05:0.001:0.15;
% Hs = 0.10:0.001:0.20;
% Hs = 0.15:0.001:0.25;
Hs = 0.2:0.001:0.25;
% Hs = 0.25:0.001:0.30;
% Hs = 0.30:0.001:0.35;
% Hs = 0.35:0.001:0.45;
data = [];

for H = Hs
    bridge = MakeBridge(R1,R2,H,theta1,theta2, V, sigma);
    res = FitParabolic(bridge);
    
    pp = res.p1;
    aa = res.alpha1;
    tt = theta1;
    ff = aa + tt - pi/2;
    rr = res.yc;
    cc = 2*rr - pp*rr^2;
    dd = 1 - pp*cc;
    % dd = cos(ff)^2 - pp * (2*sin(aa)*sin(aa+tt) - pp*sin(aa)^2);
    % dd = cos(ff)^2 - pp*cc;
    % data(end+1,1) = dd;
    data(end+1,1) = pp;
    
end

data = [Hs',data];

% data = [Hs',Ws'];

figure;
plot(Hs,data(:,2));
% plot(Hs,data(:,1), Hs,data(:,2));










