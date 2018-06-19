
clear;

sigma = 1.0;
R1 = 1.0;
% R1 = 2.0;
R2 = 1.0;

% theta1 = deg2rad(0);
% theta1 = deg2rad(20);
theta1 = deg2rad(40);
% theta2 = deg2rad(0);
% theta2 = deg2rad(20);
% theta2 = deg2rad(30);
theta2 = deg2rad(40);
% theta2 = deg2rad(60);

H = 0.2;
% H = 0.3;
% H = 0.4;
H = 0.274;

% volume large enough
V = H^3 * 10;

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
% RunEvolver;

%
[~,guess] = BridgeForceHR2(R1,R2,H,theta1,theta2,V,sigma);
alpha1 = guess.alpha1;
alpha2 = guess.alpha2;
pres = guess.pres;
vol = V;

alpha_tol = 1.0e-5;

data = [];

alpha1_hi = alpha1;
alpha1_lo = 0;



for step = 1:10000
    
    if mod(step,1) == 0
        disp(['step=',int2str(step),';alo=',num2str(alpha1_lo),';ahi=',num2str(alpha1_hi)]);
    end
    
    alpha1 = (alpha1_hi+alpha1_lo)/2;
    alpha2 = alpha1;
    pres = pres;
    
    [alpha1,alpha2,pressol, exitflag] = ParamSolve2(bridge, alpha1,alpha2,pres);
    
    if exitflag < 1
        alpha1_lo = alpha1;
    else
        alpha1_hi = alpha1;
        pres = pressol;
        
        
        %
        M = -pres / 2;
        C = R1*sin(alpha1)*sin(alpha1+theta1) + (R1^2)*M*sin(alpha1)^2;
        
        phi1 = - (pi/2 - alpha1 - theta1);
        phi2 = pi/2 - alpha2 - theta2;
        
        if mod(step,1) == 0
            phis = linspace(phi1,phi2, 50);
            [xs,ys] = ParamCurve(phis, phi1, C,M);
            % remember to shift x
            xs = xs + R1*cos(alpha1);
            
            figure(hfig);
            hold on;
            plot(xs,ys,'.-b', [xs(1),xs(end)],[ys(1),ys(end)],'xr');
            hold off;
            drawnow;
            % pause;
        end
        
        
        vol = ParamVolume(phi1,phi2,C,M);
        vol1 = SphereCapVolume(R1*sin(alpha1),R1*(1-cos(alpha1)));
        vol2 = SphereCapVolume(R2*sin(alpha2),R2*(1-cos(alpha2)));
        vol = vol-vol1-vol2;
        
        % save data
        data(end+1,:) = [alpha1,vol,alpha2,pres];
    end
    
    
    if alpha1_hi-alpha1_lo <= alpha_tol
        break;
    end
end




