
function [alpha1,alpha2,pres] = ParamSolvePres(bridge, alpha1, alpha2, pres)

R1 = bridge.R1;
R2 = bridge.R2;
theta1 = bridge.theta1;
theta2 = bridge.theta2;
H = bridge.H;
V = bridge.V;

fun = @(unk) residfunc(bridge, unk(1),unk(2),unk(3));

[sol,resid,exitflag] = fsolve(fun, [alpha1;alpha2;pres]);
disp(['exitflag=',int2str(exitflag)]);

alpha1 = sol(1);
alpha2 = sol(2);
pres = sol(3);


return
end

function [resid] = residfunc(bridge, alpha1,alpha2,pres)
    R1 = bridge.R1;
    R2 = bridge.R2;
    theta1 = bridge.theta1;
    theta2 = bridge.theta2;
    H = bridge.H;
    V = bridge.V;
    
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
    
    
    vol1 = ParamVolume(phi1,0,C1,M);
    vol2 = ParamVolume(phi2,0,C2,M);
    vcap1 = SphereCapVolume(r1,d1);
    vcap2 = SphereCapVolume(r2,d2);
    vol = vol1 + vol2 - vcap1 - vcap2;
    
    % residual
    resid = zeros(3,1);
    resid(1) = C1 - C2;
    resid(2) = (x1 + x2) - (H + d1 + d2);
    resid(3) = vol/V - 1;
    
return
end






