
function [alpha1,alpha2,pres,exitflag] = ParamSolve2(bridge, alpha1, alpha2, pres)

R1 = bridge.R1;
R2 = bridge.R2;
theta1 = bridge.theta1;
theta2 = bridge.theta2;
H = bridge.H;
% V = bridge.V;

% fun = @(unk) residfunc(bridge, alpha1, unk(1),unk(2));
fun = @(unk) residfunc2(bridge, alpha1, unk(1),unk(2));

tol = 1.0e-12;
options = optimset('TolFun',tol, 'TolX',tol, 'Display','off');
[sol,resid,exitflag] = fsolve(fun, [alpha2;pres], options);
disp(['exitflag=',int2str(exitflag)]);
resid

if norm(resid) > 1.0e-8
    exitflag = -99;
end

alpha2 = sol(1);
pres = sol(2);


return
end

function [resid] = residfunc(bridge, alpha1, alpha2,pres)
    R1 = bridge.R1;
    R2 = bridge.R2;
    theta1 = bridge.theta1;
    theta2 = bridge.theta2;
    H = bridge.H;
    
    phi1 = -(pi/2-alpha1-theta1);
    phi2 = -(pi/2-alpha2-theta2);
    
    M = -pres / 2;
    C1 = R1*sin(alpha1)*sin(alpha1+theta1) + (R1^2)*M*sin(alpha1)^2;
    C2 = R2*sin(alpha2)*sin(alpha2+theta2) + (R2^2)*M*sin(alpha2)^2;
    
    % neck points
    [x1,y1] = ParamCurve(0,phi1, C1,M);
    [x2,y2] = ParamCurve(0,phi2, C2,M);
    
    % immersed height
    r1 = R1*sin(alpha1);
    r2 = R2*sin(alpha2);
    d1 = R1*(1-cos(alpha1));
    d2 = R2*(1-cos(alpha2));
    
    
    % vol1 = ParamVolume(phi1,0,C1,M);
    % vol2 = ParamVolume(phi2,0,C2,M);
    % vcap1 = SphereCapVolume(r1,d1);
    % vcap2 = SphereCapVolume(r2,d2);
    % vol = vol1 + vol2 - vcap1 - vcap2;
    
    % residual
    resid = zeros(2,1);
    resid(1) = C1 - C2;
    resid(2) = (x1 + x2) - (H + d1 + d2);
    % resid(3) = vol/V - 1;
    
return
end

function [resid] = residfunc2(bridge, alpha1, alpha2,pres)
    R1 = bridge.R1;
    R2 = bridge.R2;
    theta1 = bridge.theta1;
    theta2 = bridge.theta2;
    H = bridge.H;
    X1 = bridge.X1;
    X2 = bridge.X2;
    
    phi1 = -(pi/2-alpha1-theta1);
    phi2 = -(pi/2-alpha2-theta2);
    
    M = -pres / 2;
    C1 = R1*sin(alpha1)*sin(alpha1+theta1) + (R1^2)*M*sin(alpha1)^2;
    C2 = R2*sin(alpha2)*sin(alpha2+theta2) + (R2^2)*M*sin(alpha2)^2;
    
    % neck points
    [x2,y2] = ParamCurve(-phi2,phi1, C1,M);
    
    % immersed height
    r1 = R1*sin(alpha1);
    r2 = R2*sin(alpha2);
    d1 = R1*(1-cos(alpha1));
    d2 = R2*(1-cos(alpha2));
    
    
    % vol1 = ParamVolume(phi1,0,C1,M);
    % vol2 = ParamVolume(phi2,0,C2,M);
    % vcap1 = SphereCapVolume(r1,d1);
    % vcap2 = SphereCapVolume(r2,d2);
    % vol = vol1 + vol2 - vcap1 - vcap2;
    
    % residual
    resid = zeros(2,1);
    resid(1) = C1 - C2;
    resid(2) = sqrt((x2+R1*cos(alpha1)-X2).^2 + y2.^2) - R2;
    % resid(3) = vol/V - 1;
    
return
end









