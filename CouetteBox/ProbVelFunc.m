
function [ u,v ] = ProbVelFunc(x,y)

EBGlobals;

r = sqrt(x.^2 + y.^2);

if r < 0.5*(R0+R1)
    u = -Omega0 * y;
    v =  Omega0 * x;
else
    u = 0;
    v = 0;
end


return
end




