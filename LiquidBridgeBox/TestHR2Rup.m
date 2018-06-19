

clear;

sigma = 1.0;
R = 1.0;
theta = deg2rad(0);
% theta = deg2rad(10);
% theta = deg2rad(20);
% theta = deg2rad(30);
% theta = deg2rad(40);
theta1 = theta;
theta2 = theta;

% V = 8.0e-5;
% V = 8.0e-4;
% V = 8.0e-3;
% V = 1.0e-2;
% V = 8.0e-2;
V = 5.0e-1;

% H = 0.0;

% bridge = MakeBridge(R1,R2,H,theta1,theta2,V,sigma);
% ff = AxisymEvolverDirectForce(bridge, 50,[],[]);

Hrup = BridgeRuptureLian(theta,V);


data = [];
num = 100;
for H = linspace(Hrup*0.8,Hrup*1.2,num)
    [F,res] = BridgeForceHR2(R,R,H,theta,theta,V,sigma);
    
    % % droplet area
    % vv = V * 0.5;
    % bb = res.b1;
    % dd = res.d1;
    % hh = fsolve(@(tt) tt.^3 + 3*bb^2*tt - dd*(3*bb^2+dd^2) - 6/pi*vv, H/2);
    % rr = (hh^2 + bb^2) / (2*hh);
    % ss = 2*pi*rr*hh;
    
    % % bridge area
    % rho1 = res.rout;
    % rho2 = res.rin;
    % phi1 = pi/2 - res.alpha1 - theta1;
    % phi2 = pi/2 - res.alpha2 - theta2;
    % ll = 2*pi*rho1 * ((rho1+rho2)*(phi1+phi2) - rho1*(sin(phi1)+sin(phi2)));
    
    % data(end+1,:) = [H, F, res.alpha1, ll-ss-ss ];
    
    data(end+1,:) = [H, F, res.alpha1 ];
end

figure;
% plot(data(:,1),data(:,2),'-', data(:,1),data(:,3),'.-');
% legend('F','alpha');
plot(data(:,1),data(:,3),'.-');
legend('alpha');
% plot(data(:,1),data(:,4),'.-');
% legend('area');
% grid on;

hold on;
plot([Hrup,Hrup],[min(data(:,3)),max(data(:,3))],'x-');
% plot([Hrup,Hrup],[min(data(:,4)),max(data(:,4))],'x-');
hold off;

[min_alpha,min_ind] = min(data(:,3));
Hmod = data(min_ind,1);
disp(['Hrup=',num2str(Hrup),';Hmod=',num2str(Hmod)]);


