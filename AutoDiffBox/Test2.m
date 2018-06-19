
clear;

% afun = @(x) tan(x) - 1 - x;
% [sol,ok] = SolveFunc(afun, 0.8, 1.0e-12);


R1 = 1.0;
R2 = 1.0;
H = 0.0;
theta1 = 0;
theta2 = 0;
sigma = 1.0;

V = 1.0e-2;
% V = 1.0e-4;
% V = 1.0e-6;
% V = 1.0e-8;
% V = 1.0e-10;


F = BridgeForceHR2(R1,R2,H,theta1,theta2,V,sigma);

