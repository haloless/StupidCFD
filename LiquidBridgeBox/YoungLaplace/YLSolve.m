function [xs,zs, res, exitflag] = YLSolve(bridge, alpha1, pres)

R1 = bridge.R1;
R2 = bridge.R2;
H = bridge.H;
theta1 = bridge.theta1;
theta2 = bridge.theta2;
V = bridge.V;


%
alpha1_lo = nan;
alpha1_hi = nan;

% 
% if isempty
% [~,guess] = BridgeForceHR2(R1,R2,H,theta1,theta2,V,sigma);
% alpha1 = guess.alpha1;
% pguess = guess.pres;

max_step = 100;

for step = 1:max_step
    
    %
    [xs,zs, res] = YLSolvePres(bridge, alpha1, pres);
    alpha2 = res.alpha2;
    
    % calc volume
    rs = zs(:,1);
    np = numel(xs);
    
    rm = 0.5 * (rs(1:end-1) + rs(2:end));
    xl = xs(2:end) - xs(1:end-1);
    vol = sum(pi .* rm.^2 .* xl);
    vol1 = SphereCapVolume(R1*sin(alpha1),R1*(1-cos(alpha1)));
    vol2 = SphereCapVolume(R2*sin(alpha2),R2*(1-cos(alpha2)));
    vol = vol - vol1 - vol2;
    
    volerr = vol - V;
    
    disp(['step=',int2str(step),';|vol|=',num2str(volerr)]);
    
    if abs(volerr)/V < 1.0e-6
        disp('volume converged');
        break;
    else
        if vol > V
            alpha1_hi = alpha1;
        else
            alpha1_lo = alpha1;
        end
        
        if isnan(alpha1_lo)
            alpha1 = alpha1 * 0.99;
        elseif isnan(alpha1_hi)
            alpha1 = alpha1 * 1.01;
        else
            alpha1 = 0.5 * (alpha1_hi + alpha1_lo);
        end
    end
end




return
end


