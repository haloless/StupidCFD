function [sol] = YLShootBVP(bridge, guess)
    R1 = bridge.R1;
    R2 = bridge.R2;
    H = bridge.H;
    theta1 = bridge.theta1;
    theta2 = bridge.theta2;
    V = bridge.V;
    
    if R2 <= 0
        R2 = 0;
    end
    
    function [dvoldx] = yl_dvol(x,y)
        ys2 = 0;
        if x < R1
            ys2 = R1^2 - x.^2;
        elseif x > R1+H
            ys2 = R2^2 - (R1+H+R2 - x).^2;
        else
            ys2 = 0;
        end
        
        dvoldx = pi * (y.^2 - ys2);
    return
    end
    
    function [dz] = yl_ode(t,z, param)
        alpha1 = param(1);
        alpha2 = param(2);
        pres = param(3);
        M = -pres / 2;
        
        ca1 = cos(alpha1);
        ca2 = cos(alpha2);
        L = R1*(1-ca1) + H + R2*(1-ca2);
        
        x = t*L + R1*ca1;
        y = z(1);
        dy = z(2);
        q = 1 + dy.^2;
        
        dz = zeros(3,1);
        % d(y)/dt
        dz(1) = L * dy;
        % d(ydot)/dt
        dz(2) = L * (q./y + 2*M*sqrt(q.^3));
        % d(v)/dt
        dz(3) = L * yl_dvol(x,y);
    return
    end
    
    function [bc] = yl_bc(z1,z2, param)
        alpha1 = param(1);
        alpha2 = param(2);
        pres = param(3);
        M = -pres / 2;
        
        ca1 = cos(alpha1);
        ca2 = cos(alpha2);
        L = R1*(1-ca1) + H + R2*(1-ca2);
        
        x1 = R1*ca1;
        x2 = R1*ca1 + L;
        
        b1 = R1 * sin(alpha1);
        phi1 = -(pi/2 - alpha1 - theta1);
        if R2 > 0
            b2 = R2 * sin(alpha2);
            phi2 = pi/2 - alpha2 - theta2;
        else
            b2 = alpha2;
            phi2 = pi/2 - theta2;
        end
        
        y1 = z1(1);
        dy1 = z1(2);
        vol1 = z1(3);
        y2 = z2(1);
        dy2 = z2(2);
        vol2 = z2(3);
        
        bc = zeros(6,1);
        % endpoints on spheres
        bc(1) = y1 - b1;
        bc(2) = y2 - b2;
        % contact angles
        % bc(3) = dy1 - tan(phi1);
        % bc(4) = dy2 - tan(phi2);
        bc(3) = atan2(dy1,1) - (phi1);
        bc(4) = atan2(dy2,1) - (phi2);
        % volume
        bc(5) = vol1 / V;
        bc(6) = vol2 / V - 1;
    return
    end
    
    function [resid, ts,zs] = yl_shoot(param)
        alpha1 = param(1);
        alpha2 = param(2);
        pres = param(3);
        M = -pres / 2;
        
        y1 = R1 * sin(alpha1);
        dy1 = -tan(pi/2 - alpha1 - theta1);
        v1 = 0;
        z1 = [ y1, dy1, v1 ];
        
        %
        odefun = @(t,z) yl_ode(t,z, param);
        % options = odeset('RelTol',1.0e-4, 'AbsTol',1.0e-6, 'BDF','on');
        options = odeset('RelTol',1.0e-9, 'AbsTol',1.0e-9);
        % options = odeset('RelTol',1.0e10, 'AbsTol',1.0e10, 'InitialStep',5e-3,'MaxStep',5e-3, ...
        % 'BDF','on', 'MaxOrder',2);
        [ts,zs] = ode45(odefun, [0,1], z1, options);
        % [ts,zs] = ode15s(odefun, [0,1], z1, options);
        
        %
        z2 = zs(end,:);
        
        % 
        bc = yl_bc(z1,z2, param);
        
        resid = zeros(3,1);
        resid(1) = bc(2);
        resid(2) = bc(4);
        resid(3) = bc(6);
        
    return
    end
    
    if ~exist('guess','var') || isempty(guess)
        if 1
            % use parabolic fitting as initial guess
            guess = FitParabolic(bridge);
            alpha1 = guess.alpha1;
            alpha2 = guess.alpha2;
            pres = guess.p1;
        else
            % use circular approx.
            [~,guess] = BridgeForceHR2(R1,R2,H,theta1,theta2,V, 1.0);
            alpha1 = guess.alpha1;
            if R2 > 0
                alpha2 = guess.alpha2;
            else
                alpha2 = guess.b2;
            end
            pres = guess.pres;
        end
        % [alpha1,alpha2,pres]
    else
        alpha1 = guess(1);
        alpha2 = guess(2);
        pres = guess(3);
    end
    
    %
    options = optimset('TolFun',1.0e-6, 'TolX',1.0e-8);
    % [par,resid,exitflag] = fsolve(@yl_shoot, [alpha1;alpha2;pres], options);
    [par,resid,exitflag] = mysolve(@yl_shoot, [alpha1;alpha2;pres], [1.0e-6;1.0e-5;1.0e-5]);
    
    %
    sol = struct();
    sol.flag = exitflag;
    sol.resid = resid;
    sol.alpha1 = par(1);
    sol.alpha2 = par(2);
    sol.pres = par(3);
    
    %
    [~,ts,zs] = yl_shoot(par);
    if R2 > 0
        L = R1*(1-cos(sol.alpha1)) + H + R2*(1-cos(sol.alpha2));
    else
        L = R1*(1-cos(sol.alpha1)) + H;
    end
    sol.xs = R1*cos(sol.alpha1) + L.*ts;
    sol.ys = zs(:,1);
return
end



function [x,f,exitflag] = mysolve(fun, xguess, tol)
    
    n = numel(xguess);
    
    h = 1.0e-6;
    
    x = xguess;
    f = fun(x);
    
    exitflag = 0;
    
    for iter = 1:200
        
        df = zeros(n,n);
        for i = 1:n
            xp = x;
            xp(i) = xp(i) + h;
            fp = fun(xp);
            
            xm = x;
            xm(i) = xm(i) - h;
            fm = fun(xm);
            
            df(:,i) = (fp-fm) ./ (h*2);
        end
        
        dx = df \ f;
        
        x = x - dx * 0.4;
        f = fun(x);
        
        if all(abs(f)<=tol)
            exitflag = 1;
            break;
        end
    end
    
    return
end








