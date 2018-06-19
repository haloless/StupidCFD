function [sol] = YLSolveBVP(bridge)
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
        bc(3) = dy1 - tan(phi1);
        bc(4) = dy2 - tan(phi2);
        % volume
        bc(5) = vol1 / V;
        bc(6) = vol2 / V - 1;
    return
    end
    
    % use parabolic fitting as initial guess
    guess = FitParabolic(bridge);
    alpha1 = guess.alpha1;
    alpha2 = guess.alpha2;
    pres = guess.p1;
    % [alpha1,alpha2,pres]
    
    function [zinit] = yl_init(t)
        ca1 = cos(alpha1);
        ca2 = cos(alpha2);
        L = R1*(1-ca1) + H + R2*(1-ca2);
        
        a = guess.a;
        b = guess.b;
        c = guess.c;
        
        x = R1*ca1 + L*t;
        y = a*x^2 + b*x + c;
        dy = 2*a*x + b;
        vol = t * V;
        
        zinit = zeros(3,1);
        zinit(1) = y;
        zinit(2) = dy;
        zinit(3) = vol;
    return
    end
    
    % with initial mesh
    solinit = bvpinit(linspace(0, 1, 201), @yl_init, [alpha1, alpha2, pres]);
    
    % options
    options = bvpset('RelTol',1.0e-9, 'AbsTol',1.0e-9, 'Stats','off');
    
    % solution 
    sol = bvp4c(@yl_ode, @yl_bc, solinit, options);
    
    % check convergence
    % this is a hack: BVP4C seems to return a converged solution with a STATS field
    if isfield(sol, 'stats')
        sol.ok = 1;
    else
        sol.ok = 0;
    end
    
return
end







