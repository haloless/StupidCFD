
function TestAxisymEvolverRup(V, H0, theta1,theta2)
% Find the rupture distance empirically

% clear all;

sigma = 1.0;
R1 = 1.0;
% R2 = 1.0;
% R2 = 2.0;
R2 = -1.0;

% R1 = 2.381;
% R2 = 2.381;
% R2 = 1.588;
% R2 = 1.191;
% R2 = -1;

Rm = DerjaguinRadius(R1,R2);

if ~exist('theta1','var')
    % theta1 = deg2rad(0);
    % theta1 = deg2rad(10);
    % theta1 = deg2rad(15);
    % theta1 = deg2rad(20);
    theta1 = deg2rad(30);
    % theta1 = deg2rad(40);
    % theta1 = deg2rad(60);
else
    theta1 = deg2rad(theta1);
end
if ~exist('theta2','var')
    % theta2 = deg2rad(0);
    % theta2 = deg2rad(10);
    % theta2 = deg2rad(15);
    % theta2 = deg2rad(20);
    % theta2 = deg2rad(30);
    % theta2 = deg2rad(40);
    % theta2 = deg2rad(50);
    theta2 = deg2rad(60);
else
    theta2 = deg2rad(theta2);
end


% V = 1.0e-5;
% V = 1.0e-4;
% V = 1.0e-4;
% V = 1.0e-3;
% V = 1.0e-2;
% V = 3.0e-1;

if isempty(H0)
    % H0 = BridgeRuptureLian((theta1+theta2)/2, V);
    H0 = BridgeRuptureLian(min(theta1,theta2), V/Rm^3) * Rm;
    % H0 = H0 *0.8;
    H0 = floor(H0/0.01) * 0.01;
end


hfig = figure;


H = H0;
Hincr = 0.001;

while 1
    
    bridge = MakeBridge(R1,R2,H,theta1,theta2,V,sigma);

    %
    X1 = bridge.X1;
    X2 = bridge.X2;
    
    clf(hfig);
    
    % draw spheres
    hold on;
    PlotBridgeGeom(bridge);
    hold off;
    
    % use odd number 
    np = 51;
    % np = 81;
    % np = 101;
    
    % initialize
    [rp,xp] = AxisymEvolverGuessInit(bridge,np);
    if 0
        % draw initial shape
        figure(hfig);
        hold on;
        plot(xp,rp,'-b');
        hold off;
        drawnow;
    end
    
    % optimize
    [rp,xp,pres,sene] = AxisymEvolver(bridge,np,rp,xp);

    
    
    if 1
        % log current solution
        % force
        [F1,~,~,F2,~,~] = AxisymEvolverEvalForce(bridge, np,rp,xp,pres)
        % energy
        sene;
        alpha1 = asin(rp(1)/R1)
        alpha2 = asin(rp(end)/R2);
        pres
        plim1 = sin(alpha1+theta1)/sin(alpha1);
        plim2 = sin(alpha2+theta2)/sin(alpha2);
    end
    
    if 1
        % volume distribution
        % always calculate the left side, i.e. sphere 1
        
        % find the neck index
        [~,ineck] = min(rp(:));
        ii = 1:(ineck-1);
        
        xneck = xp(ineck)
        (xneck-R1)/H
        rneck = rp(ineck)
        
        vfrac1 = sum(pi * (0.5*(rp(ii)+rp(ii+1))).^2 .* (xp(ii+1)-xp(ii)));
        vcap1 = SphereCapVolume(R1*sin(alpha1),R1*(1-cos(alpha1)));
        vfrac1 = (vfrac1-vcap1) / V
        
        ii = ineck:(np-1);
        vfrac2 = sum(pi * (0.5*(rp(ii)+rp(ii+1))).^2 .* (xp(ii+1)-xp(ii)));
        vfrac2 = vfrac2 / V;
        
        b2 = rp(end)
        
        % draw neck position
        figure(hfig);
        hold on;
        plot([xp(ineck),xp(ineck)],[0,rp(ineck)],'x-');
        hold off;
    end
    
    if 1
        % draw optimized shape
        figure(hfig);
        hold on;
        plot(xp,rp,'.-r');
        hold off;
        title(['V=',num2str(V), '; f1=',num2str(vfrac1), '; H=',num2str(H),'; Hold=',num2str(H-Hincr)]);
        drawnow;
    end 
    
    
    % csvwrite('tmp.csv', [xp,rp]);
    
    hh = input(['incr(',num2str(Hincr),')=']);
    if ~isempty(hh)
        Hincr = hh;
    end
    H = H + Hincr;
    
    % data = [xp,rp];
    
    
end



return
end







