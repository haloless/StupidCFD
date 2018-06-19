
function [ res ] = FitParabolic(bridge)
    
    % bridge parameters
    R1 = bridge.R1;
    R2 = bridge.R2;
    theta1 = bridge.theta1;
    theta2 = bridge.theta2;
    H = bridge.H;
    V = bridge.V;
    
    if R2 <= 0
        R2 = 0;
    end
    
    % 
    func = @(x) residfunc(bridge, x(1),x(2),x(3),x(4),x(5));
    
    % guess = [ 0; 0; R1/2; pi/6; pi/6 ];
    [bguess, a1guess,a2guess] = SolveStraightBridge(R1,R2,H,V);
    if R2 > 0
        guess = [ 0; 0; bguess; a1guess; a2guess ];
    else
        guess = [ 0; 0; bguess; a1guess; bguess ];
    end
    
    options = optimoptions('fsolve','TolFun',1.0e-4,'TolX',1.0e-5, ...
    'MaxFunEvals',5000, 'MaxIter',10000);
    [sol,val,exitflag] = fsolve(func, guess, options);
    
    if exitflag ~= 1
        warning('Failed to solve parabolic');
        % failed? try return straight guess instead
        sol = guess;
    end
    
    
    % more info
    a = sol(1);
    b = sol(2);
    c = sol(3);
    alpha1 = sol(4);
    alpha2 = sol(5);
    
    ca1 = cos(alpha1);
    sa1 = sin(alpha1);
    ca2 = cos(alpha2);
    sa2 = sin(alpha2);
    
    % position of contact points
    x1 = R1 * ca1;
    x2 = R1 + H + R2 - R2 * ca2;
    [y1,dy1] = parafunc(a,b,c, x1);
    [y2,dy2] = parafunc(a,b,c, x2);
    
    % immersed height
    h1 = R1 * (1-ca1);
    h2 = R2 * (1-ca2);
    
    vcap1 = pi/3 * R1^3 * (2 - 3*ca1 + ca1.^3);
    vcap2 = pi/3 * R2^3 * (2 - 3*ca2 + ca2.^3);
    
    scap1 = 2*pi*R1 * h1;
    scap2 = 2*pi*R2 * h2;
    
    % area
    spara = integral(@(x) surffunc(a,b,c,x), x1,x2);
    sene = spara - cos(theta1)*scap1 - cos(theta2)*scap2;
    
    %
    p1 = presfunc(a,b,c, x1);
    p2 = presfunc(a,b,c, x2);
    xc = -b ./ (2*a);
    yc = parafunc(a,b,c, xc);
    pc = presfunc(a,b,c, xc);
    
    % F1 = 2*pi*R1*sin(alpha1)*sin(alpha1+theta1) - p1*pi*R1^2*sin(alpha1)^2;
    % F2 = 2*pi*R2*sin(alpha2)*sin(alpha2+theta2) - p2*pi*R2^2*sin(alpha2)^2;
    % Fc = 2*pi*yc - (p1+p2)/2 *pi*yc^2;
    F1 = 2*pi*R1*sin(alpha1)*sin(alpha1+theta1);
    F2 = 2*pi*R2*sin(alpha2)*sin(alpha2+theta2);
    Fc = 2*pi*yc;
    
    % save final results
    res = struct();
    % solver state
    res.resid = val;
    res.exitflag = exitflag;
    % solution
    res.a = a;
    res.b = b;
    res.c = c;
    res.alpha1 = alpha1;
    res.alpha2 = alpha2;
    % 
    res.x1 = x1;
    res.y1 = y1;
    res.dy1= dy1;
    res.x2 = x2;
    res.y2 = y2;
    res.dy2 = dy2;
    % 
    res.h1 = h1;
    res.h2 = h2;
    res.scap1 = scap1;
    res.scap2 = scap2;
    res.spara = spara;
    %
    res.W = sene;
    res.WL = spara;
    res.W1 = -cos(theta1)*scap1;
    res.W2 = -cos(theta2)*scap2;
    
    %
    res.p1 = p1;
    res.p2 = p2;
    res.pc = pc;
    res.xc = xc;
    res.yc = yc;
    res.F1 = F1;
    res.F2 = F2;
    res.Fc = Fc;
    
    
    
return
end


function [ resid ] = residfunc(bridge, a,b,c,alpha1,alpha2)
    
    R1 = bridge.R1;
    R2 = bridge.R2;
    theta1 = bridge.theta1;
    theta2 = bridge.theta2;
    H = bridge.H;
    V = bridge.V;
    
    if R2 <= 0
        R2 = 0;
    end
    
    ca1 = cos(alpha1);
    sa1 = sin(alpha1);
    ca2 = cos(alpha2);
    sa2 = sin(alpha2);
    
    x1 = R1 * ca1;
    f1 = R1 * sa1;
    g1 = tan(alpha1 + theta1 - pi/2);
    x2 = R1 + H + R2 - R2 * ca2;
    if R2 > 0
        f2 = R2 * sa2;
        g2 = tan(pi/2 - alpha2 - theta2);
    else
        f2 = alpha2;
        g2 = tan(pi/2 - theta2);
    end
    
    % sphere cap
    vcap1 = pi/3 * R1^3 * (2 - 3*ca1 + ca1.^3);
    vcap2 = pi/3 * R2^3 * (2 - 3*ca2 + ca2.^3);
    
    %
    vol = volintfunc(a,b,c, x2) - volintfunc(a,b,c, x1);
    vol = vol - vcap1 - vcap2;
    % volumetric error
    volerr = vol./V - 1;
    % volerr = (vol./V).^(1/3) - 1;
    
    %
    [y1,dy1] = parafunc(a,b,c, x1);
    [y2,dy2] = parafunc(a,b,c, x2);
    
    % residual
    resid = [ y1-f1; y2-f2; dy1-g1; dy2-g2; volerr ];
return
end

function [ surf ] = surffunc(a,b,c, x)
    [y,dy] = parafunc(a,b,c, x);
    
    surf = 2*pi .* y .* sqrt(1+dy.^2);
return
end


function [ y,dy ] = parafunc(a,b,c, x)
    y = a .* x.^2 + b .* x + c;
    dy = 2*a .* x + b;
return
end


function [ vol ] = volintfunc(a,b,c, x)
    x5 = x.^5;
    x4 = x.^4;
    x3 = x.^3;
    x2 = x.^2;
    vol = (a.^2)/5 .* x5 + (a.*b)/2 .* x4 + (b.^2+2*a.*c)/3 .* x3 + (b.*c) .* x2 + (c.^2) .* x;
    vol = vol * pi;
return
end

function [ pres ] = presfunc(a,b,c, x)
    ddy = 2 * a;
    [y,dy] = parafunc(a,b,c, x);
    q = sqrt(1 + dy.^2);
    pres = 1 ./ (y.*q) - ddy./(q.^3);
return
end







