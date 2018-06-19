function [res] = BridgeIntegFindAlfa(bridge, pres, a1guess,a2guess,Hguess)
    
    R1 = bridge.R1;
    R2 = bridge.R2;
    % H = bridge.H;
    theta1 = bridge.theta1;
    theta2 = bridge.theta2;
    V = bridge.V;
    % X1 = bridge.X1;
    % X2 = bridge.X2;
    if R2 <= 0
        R2 = 0;
    end
    
    % result struct
    res = struct;

    
    % ODE function
    function [dz] = yl_ode(t,z, par)
        alpha1 = par(1);
        alpha2 = par(2);
        H = par(3);
        
        M = -pres / 2;
        
        ca1 = cos(alpha1);
        ca2 = cos(alpha2);
        sa1 = sin(alpha1);
        sa2 = sin(alpha2);
        
        L = R1*(1-ca1) + H + R2*(1-ca2);
        x = R1*ca1 + L*t;
        
        y = z(1);
        dy = z(2);
        q = sqrt(1+dy.^2);
        ddy = q.^2 ./ y + 2*M * q.^3;
        
        %
        ys2 = 0;
        if x < R1
            ys2 = R1^2 - x^2;
        elseif x > R1+H
            ys2 = R2^2 - (x-R1-H-R2)^2;
        end
        dvol = pi * (y^2 - ys2);
        
        dz = zeros(3,1);
        dz(1) = L * dy;
        dz(2) = L * ddy;
        dz(3) = L * dvol;
    return
    end
    
    function [resid, ts,zs, te,ze] = yl_shoot(par, varargin)
        alpha1 = par(1);
        alpha2 = par(2);
        H = par(3);
        
        y1 = R1 * sin(alpha1);
        dy1 = -tan(pi/2-alpha1-theta1);
        vol1 = 0;
        
        options = odeset('RelTol',1.0e-9, 'AbsTol',1.0e-9, varargin{:});
        % odesol = ode45(@(t,z) yl_ode(t,z,par), [0,1], [y1,dy1,vol1], options);
        % ts = odesol.x;
        % zs = odesol.y;
        if length(varargin) > 0
            [ts,zs, te,ze] = ode45(@(t,z) yl_ode(t,z,par), [0,1], [y1,dy1,vol1], options);
        else
            [ts,zs] = ode45(@(t,z) yl_ode(t,z,par), [0,1], [y1,dy1,vol1], options);
            te = [];
            ze = [];
        end
        
        
        y2 = zs(end,1);
        dy2 = zs(end,2);
        vol2 = zs(end,3);
        
        resid = zeros(3,1);
        if R2 > 0
            resid(1) = y2 - R2*sin(alpha2);
            resid(2) = dy2 - tan(pi/2-alpha2-theta2);
        else
            resid(1) = y2 - alpha2;
            resid(2) = dy2 - tan(pi/2-theta2);
        end
        resid(3) = vol2/V - 1;
    return
    end
    

    
    
    % solve
    options = optimset('TolFun',1.0e-6, 'TolX',1.0e-8, 'Display','off');
    [sol,resid,exitflag] = fsolve(@yl_shoot, [a1guess; a2guess; Hguess], options);
    % disp(['flag=',int2str(exitflag)]);
    % disp(num2str(resid));
    
    % the final solution
    alpha1 = sol(1);
    alpha2 = sol(2);
    H = sol(3);
    L = R1*(1-cos(alpha1)) + H + R2*(1-cos(alpha2));
    
    %
    function [value,isterminal,direction] = yl_event_neck(t,z)
        y = z(1);
        dy = z(2);
        vol = z(3);
        
        M = -pres/2;
        q = sqrt(1+dy.^2);
        ddy = q.^2 ./ y + 2*M * q.^3;
        
        value = zeros(3,1);
        isterminal = zeros(3,1);
        direction = zeros(3,1);
        
        %
        value(1) = dy;
        value(2) = ddy;
    return
    end
    
    % compute final bridge 
    % [resid, ylsol] = yl_shoot([alpha1; alpha2; H], 'Events',@yl_event_neck);
    % ts = ylsol.x;
    % zs = ylsol.y;
    [resid, ts,zs, te,ze] = yl_shoot([alpha1; alpha2; H], 'Events',@yl_event_neck);
    ylsol = [];
    
    %
    res.resid = resid;
    res.ts = ts;
    res.xs = ts.*L + R1*cos(alpha1);
    res.ys = zs(:,1);
    %
    res.alpha1 = alpha1;
    res.alpha2 = alpha2;
    res.H = H;
    res.pres = pres;
    res.b1 = R1*sin(alpha1);
    if R2 > 0
        res.b2 = R2 * sin(alpha2);
    else
        res.b2 = alpha2;
    end
    
    % locate neck
    [~,res.ineck] = min(res.ys);
    res.xneck = res.xs(res.ineck);
    res.yneck = res.ys(res.ineck);
    
    % events
    if ~isempty(te)
        res.xe = te .* L + R1*cos(alpha1);
        res.ye = ze(:,1);
    else
        res.xe = [];
        res.ye = [];
    end
    
return
end


