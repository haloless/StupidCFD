
function [ cmn ] = cm_stokesd(R,sgn,np,rhs)


A = zeros(np,np);

for k = 1:np
    n = k;
    fn1 = (n+1) / (4*n-2);
    fn2 = (n-2) / (n*(4*n-2));
    
    for i = 1:np
        theta = (i-1) * pi/(np-1);
        co1 = cos(theta);
        si1 = sqrt(1.0-co1^2);
        r2 = sqrt(si1^2 + (R+co1)^2);
        co2 = -(R+co1)/r2;
        
        ppns1 = cm_pns(n,1,co1);
        ppns2 = cm_pns(n,1,co2);
        
        rr = 1.0/r2;
        A(i,k) = -ppns1 - rr*sgn*ppns2/r2^(n+1);
    end
end

sol = A \ rhs;
cmn = sol;

return
end

