function [] = TestParamRup(V, H)

% clear;

sigma = 1.0;
R1 = 1.0;
% R2 = 1.0;
R2 = -1.0;

% Rm = DerjaguinRadius(R1,R2);
Rm = R1;

% theta1 = deg2rad(0);
% theta1 = deg2rad(20);
theta1 = deg2rad(30);
% theta1 = deg2rad(40);
% theta1 = deg2rad(60);
% theta2 = deg2rad(0);
% theta2 = deg2rad(2);
% theta2 = deg2rad(10);
% theta2 = deg2rad(20);
% theta2 = deg2rad(30);
% theta2 = deg2rad(40);
theta2 = deg2rad(50);
% theta2 = deg2rad(60);

% V = 0.0001;
% V = 0.0003;
% V = 0.0006;
% V = 0.001;
% V = 0.01;
% V = 0.04;
% V = 0.06;
% V = 0.1;

% V = 0.15^3;
% H = 0.099;

V = V * Rm^3;
H = H * Rm;


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
    pres = pres;
end

if 1
    
    
    % use optimal solution
    guess = [alpha1,alpha2,pres];
    if R2 <= 0
        guess(2) = rp(end);
    end
    sol = YLShootBVP(bridge, guess);
    
    figure(hfig);
    hold on;
    plot(sol.xs,sol.ys,'k-');
    hold off;
    pause;
    
    assert(sol.flag >= 1);
    
    alpha1 = sol.alpha1;
    alpha2 = sol.alpha2;
    pres = sol.pres;
    
    if 0
        Hincr = 1.0e-4;
        while (Hincr >= 1.0e-5)
            bridgeold = bridge;
            solold = sol;
            Hold = H;
            
            H = H + Hincr;
            
            disp(['H=',sprintf('%0.8f',H)]);
            
            bridge = MakeBridge(R1,R2,H,theta1,theta2,V,sigma);
            
            % use previous good solution
            guess = [ alpha1, alpha2, pres];
            sol = YLShootBVP(bridge, guess);
            
            % if sol.flag == 1
            if sol.flag <= 3 && sol.flag>=1
                % update
                alpha1 = sol.alpha1;
                alpha2 = sol.alpha2;
                pres = sol.pres;
            else
                % recover
                disp('Failed, recover and reduce');
                H = Hold;
                sol = solold;
                bridge = bridgeold;
                Hincr = Hincr * 0.1; 
                % break;
            end
        end
    end
    
    %
    % sol.resid
    
    [~,ineck] = min(sol.ys);
    xneck = sol.xs(ineck);
    yneck = sol.ys(ineck);
    vol1 = trapz(sol.xs(1:ineck),pi*sol.ys(1:ineck).^2);
    vcap1 = SphereCapVolume(R1*sin(sol.alpha1), R1*(1-cos(sol.alpha1)));
    vol1 = vol1 - vcap1;
    
    if 1
        figure(hfig);
        hold on;
        PlotBridgeGeom(bridge);
        plot(sol.xs,sol.ys,'b-', [xneck,xneck],[0,yneck],'mx-');
        hold off;
        msgstr = sprintf('H=%0.6f;H/Rm=%0.6f;V1/V=%0.4f',H,H/Rm,vol1/V);
        title(msgstr);
        disp(msgstr);
    end
    
    if 1
        % evaluate 
        vol1a = vol1 + vcap1;
        vol2 = V - vol1;
        vol1b = vol1a - vol2;
        volimg = V + vol1b - vcap1
        % hrupimg = (1+0.5*theta1)*(volimg^(1/3)+0.1*volimg^(2/3))
        hrupimg = BridgeRuptureWillet(R1,R1, theta1, volimg);
        disp(sprintf('rupture using image = %0.4f', hrupimg/2 + H-(xneck-R1)));
    end
    
    if (theta2 < theta1+alpha1)
        ipos = 0;
        apos = 0;
        emin = 99999;
        for ip = ineck:(np-1)
            slope = (rp(ip+1)-rp(ip-1)) / (xp(ip+1)-xp(ip-1));
            % rad2deg(pi/2-atan(slope))
            ee = abs((slope)-tan(pi/2-theta1));
            if ee < emin
                apos = pi/2-atan(slope);
                ipos = ip;
                emin = ee;
            end
        end
        % rad2deg(apos)
        figure(hfig);
        hold on;
        plot([xp(ipos),xp(ipos)],[0,rp(ipos)], 'cx-');
        hold off;
        
        voltrunc = trapz(xp(1:ipos), pi*rp(1:ipos).^2);
        voltrunc = voltrunc - vcap1
        hruptrunc = (1+0.3*theta1) * (voltrunc^(1/3) - 0.23*voltrunc^(2/3));
        disp(sprintf('rupture using trunc = %0.4f', hruptrunc + H+R1-xp(ipos)));
    end
    
    if (theta2 > theta1) && false
        a2hat = theta2 - theta1;
        R2hat = rp(end) / sin(a2hat);
        vcap2 = SphereCapVolume(R2hat*sin(a2hat), R2hat*(1-cos(a2hat)));
        volhat = V - vcap2;
        hruphat = BridgeRuptureWillet(R1,R2hat,theta1,volhat);
        disp(sprintf('rupture using conv = %0.4f', hruphat));
    end
end


return
end















