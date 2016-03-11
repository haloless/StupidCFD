% Calculate derivative of PNS(n,m,x)

function [ pnsp ] = cm_pnsp(n,m,x)

n1 = n - m;

if (m == 0)
    p0 = 1.0;
    pp0 = 0.0;
    if (n > 0)
        p1 = x;
        pp1 = 1.0 - x*x;
        if (n > 1)
            pmin = p0;
            p = p1;
            pp = pp1;
            for i = 2:n1
                pplus = 2.0*x*p - pmin - (x*p-pmin)/i;
                pppls = i*pp1*p + x*pp;
                pmin = p;
                p = pplus;
                pp = pppls;
            end
            pnsp = pppls;
        else
            pnsp = pp1;
        end
    else
        pnsp = 0.0;
    end
else
    coef = 1.0;
    for i = 1:m
        coef = coef * (2*i-1);
    end
    p0 = coef;
    pp0 = -m*x*p0;
    if (n > m)
        p1 = p0*x*(2*m+1);
        pp1 = (2*m+1) * (p0*(1.0-x^2) + x*pp0);
        if (n > m+1)
            p = p1;
            pp = pp1;
            pmin = p0;
            for i = 2:n1
                pplus = 2.0*x*p - pmin + (2*m-1)*(x*p-pmin)/i;
                pppls = ((i+2*m)*(1.0-x^2) + m*x^2)*p - m*x*pplus + x*pp;
                pmin = p;
                p = pplus;
                pp = pppls;
            end
            pnsp = pppls;
        else
            pnsp = pp1;
        end
    else
        pnsp = pp0;
    end
end

return
end






