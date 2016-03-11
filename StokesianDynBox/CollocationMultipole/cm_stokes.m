
function [ amn,bmn,cmn ] = cm_stokes(R,nstar,m,sgn,np,rhs)


A = zeros(np*3,np*3);

for k = 1:np
    k2 = k + np;
    k3 = k2 + np;
    
    n = k + nstar - 1;
    fn1 = (n+1) / (4*n-2);
    fn2 = (n-2) / (n*(4*n-2));
    
    for i = 1:np
        i2 = i + np;
        i3 = i2 + np;
        
        theta = (i-1) * pi/(np-1);
        co1 = cos(theta);
        si1 = sqrt(1.0-co1*co1);
        r2 = sqrt(si1^2 + (R+co1)^2);
        co2 = -(R+co1) / r2;
        
        pns1 = cm_pns(n, m, co1);
        pns2 = cm_pns(n, m, co2);
        ppns1 = cm_pns(n, m+1, co1);
        ppns2 = cm_pns(n, m+1, co2);
        pnsp1 = cm_pnsp(n, m, co1);
        pnsp2 = cm_pnsp(n, m, co2);
        
        rz = (1.0/r2)^m;
        rr = rz / r2;
        rphi = rz * r2;
        
        A(i,k) = fn1*(co1*pns1-rz*sgn*co2*pns2/r2^n) - fn2*(pnsp1-rz*sgn*pnsp2/r2^n);
        A(i,k2) = -(n-m+1) * (cm_pns(n+1,m,co1) - rz*sgn*cm_pns(n+1,m,co2)/r2^(n+2));
        A(i,k3) = m * (pns1 - rz*sgn*pns2/r2^(n+1));
        
        A(i2,k) = (fn1+fn2*m)*(pns1+rr*sgn*pns2/r2^n) + fn2*(co1*ppns1+rr*sgn*co2*ppns2/r2^n);
        A(i2,k2) = -(n+m+1)*(pns1+rr*sgn*pns2/r2^(n+2)) - co1*ppns1 - rr*sgn*co2*ppns2/r2^(n+2);
        A(i2,k3) = -ppns1 - rr*sgn*ppns2/r2^(n+1);
        
        A(i3,k) = -m*fn2 * (pns1 + sgn*rphi*pns2/r2^n);
        A(i3,k2) = m * (pns1 + rphi*sgn*pns2/r2^(n+2));
        A(i3,k3) = pnsp1 + rphi*sgn*pnsp2/r2^(n+1);
    end
end

ndim = 3 * np;
if (m == 0) 
    ndim = 2*np;
end

sol = zeros(np*3,1);
if (m == 0)
    Idim = 1:ndim;
    sol(Idim) = A(Idim,Idim) \ rhs(Idim);
else
    sol = A \ rhs;
end

amn = zeros(np,1);
bmn = zeros(np,1);
cmn = zeros(np,1);
for i = 1:np
    i2 = i + np;
    i3 = i2 + np;
    amn(i) = sol(i);
    bmn(i) = sol(i2);
    cmn(i) = sol(i3);
end




return
end

