function [ts,zs, res] = YLIntegrateBridge(bridge, alpha1,pres)
    
    R1 = bridge.R1;
    R2 = bridge.R2;
    H = bridge.H;
    theta1 = bridge.theta1;
    theta2 = bridge.theta2;
    
    if R2 <= 0
        R2 = 0;
    end
    
    %
    x1 = R1 * cos(alpha1);
    y1 = R1 * sin(alpha1);
    dy1 = -cot(alpha1+theta1);
    
    % 
    if R2 > 0
        [alpha2,b2] = solve_alpha2(R1,R2,theta1,theta2, pres, alpha1);
        
        x2 = R1 + H + R2 - R2*cos(alpha2);
        y2 = R2 * sin(alpha2);
        dy2 = cot(alpha2+theta2);
    else
        alpha2 = 0;
        x2 = R1 + H;
        y2 = 0;
        dy2 = cot(theta2);
    end
    
    
    %
    dzfunc = @(t,z) [ z(2); ddy_func(z(1),z(2),pres)];
    options = odeset('RelTol',1.0e-8, 'AbsTol',1.0e-8, 'MaxStep',R1*1e-2);
    
    % perform integration
    [ts,zs] = ode45(dzfunc, [x1,x2], [y1,dy1]);
    % [ts,zs] = ode45(dzfunc, [x1,x2], [y1,dy1], options);
    % [ts,zs] = ode113(dzfunc, [x0,xe], [y0,dy0], options);
    
    % 
    res = struct();
    res.alpha1 = alpha1;
    res.x1 = x1;
    res.y1 = y1;
    res.dy1 = dy1;
    res.alpha2 = alpha2;
    res.x2 = x2;
    res.y2 = y2;
    res.dy2 = dy2;
return
end


function [alpha2,b2] = solve_alpha2(R1,R2,theta1,theta2, pres, alpha1)
    
    b1 = R1 * sin(alpha1);
    
    ffunc = @(rr,aa,tt,pp) 2*rr.*sin(aa).*sin(aa+tt) - rr^2*pp.*(sin(aa).^2);
    
    alpha2func = @(alpha2) ffunc(R1,alpha1,theta1,pres) - ffunc(R2,alpha2,theta2,pres);
    
    options = optimset('TolFun',1.0e-10, 'TolX',1.0e-10, 'Display','off');
    [alpha2,~,exitflag] = fsolve(alpha2func, alpha1, options);
    if (exitflag ~= 1)
        error('Failed to solve alpha2, FLAG=%d', exitflag);
    end
    
    b2 = R2 * sin(alpha2);
    
return
end

function [q] = q_func(dy)
    q = 1 + dy.^2;
return
end

function [ddy] = ddy_func(y,dy,pp)
    q = 1 + dy.^2;
    ddy = q./y - pp.*sqrt(q.^3);
return
end






