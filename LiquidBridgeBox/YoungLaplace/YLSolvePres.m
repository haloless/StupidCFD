function [ts,zs, res] = YLSolvePres(bridge, alpha1, pguess)

R1 = bridge.R1;
R2 = bridge.R2;
H = bridge.H;
theta1 = bridge.theta1;
theta2 = bridge.theta2;

% guess if not provided
if isempty(pguess)
    rhoguess = 0.5 * (R1+H+R2 - R1*cos(alpha1) - R2*cos(alpha1));
    pguess = 1/(R1*sin(alpha1)) - 1/rhoguess;
end

pres_lo = nan;
pres_hi = nan;
err_lo = 0;
err_hi = 0;

% pres = (pres_lo+pres_hi)*0.5;


pres = pguess;


% maxstep = 1;
maxstep = 100;

for step = 1:maxstep
    
    % 
    [ts,zs, res] = YLIntegrateBridge(bridge, alpha1,pres);
    
    % numerical result of end point
    x2 = ts(end);
    y2 = zs(end,1);
    dy2 = zs(end,2);
    
    tol = 1.0e-6;
    x2err = (x2 - res.x2);
    y2err = (y2 - res.y2);
    dy2err = (dy2 - res.dy2);
    
    if 0
        disp(['step=',int2str(step), ...
        ';|x2|=',num2str(x2err),';|y2|=',num2str(y2err),';|dy2|=',num2str(dy2err)]);
    end
    
    if abs(x2err)<=tol
        if abs(y2err)<=tol && abs(dy2err)<=tol
            disp('pressure converged');
            break;
        else
            if y2err < 0
                pres_hi = pres;
            else
                pres_lo = pres;
            end
            
            if isnan(pres_lo)
                % pres = pres * 1.01;
                % pres = pres + sign(pres)*0.1;
                pres = pres - 0.1;
            elseif isnan(pres_hi)
                % pres = pres * 0.99;
                % pres = pres - sign(pres)*0.1;
                pres = pres + 0.1;
            else
                pres = (pres_lo+pres_hi)*0.5;
            end
        end
    else
        % error('Failed to integrate');
        pres = pres * 0.95;
        % if pres > 0
            % pres = pres - 0.1;
        % else
            % pres = pres + 0.1;
        % end
    end
end

res.pres = pres;
% res.alpha2 = alpha2;


return
end





