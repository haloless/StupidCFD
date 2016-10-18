
function [ a ] = FmmMoment(y,node, ielem,ioff, num, nexp, cx,cy, u,uoff, bc,dnorm)

a = zeros(nexp+1,1);

for i = 1:num
    nelm = ielem(i+ioff-1);
    n1 = node(1,nelm);
    n2 = node(2,nelm);
    z1 = complex(y(1,n1)-cx, y(2,n1)-cy);
    z2 = complex(y(1,n2)-cx, y(2,n2)-cy);
    zwbar = conj(z2-z1);
    zwbar = zwbar / abs(zwbar);
    zp1 = z1 * zwbar;
    zp2 = z2 * zwbar;
    znorm = complex(dnorm(1,nelm),dnorm(2,nelm));
    
    if bc(1,nelm) == 1
        phi = 0.0;
        q = u(i+uoff-1);
    elseif bc(1,nelm) == 2
        phi = u(i+uoff-1);
        q = 0.0;
    end
    
    a(1) = a(1) - (zp2-zp1)*q;
    for k = 1:nexp
        a(k+1) = a(k+1) + (zp2-zp1)*znorm*phi;
        zp1 = zp1 * z1 / (k+1);
        zp2 = zp2 * z2 / (k+1);
        a(k+1) = a(k+1) - (zp2-zp1)*q;
    end
end

return
end


