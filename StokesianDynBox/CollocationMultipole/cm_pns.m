% Calculates P_m^n(x) [sin(theta)]^(-m)

function [ pns ] = cm_pns(n,m,x)

if (m <= n)
    n1 = n - m;
    if (m == 0) 
        p0 = 1.0;
        if (n > m)
            p1 = x;
            if (n > m+1)
                pminus = p0;
                p = p1;
                for i = 2:n1
                    pplus = 2.0*x*p - pminus - (x*p-pminus)/i;
                    pminus = p;
                    p = pplus;
                end
                pns = pplus;
            else
                pns = p1;
            end
        else
            pns = p0;
        end
    else
        coef = 1.0;
        for i = 1:m
            coef = coef * (2*i-1);
        end
        p0 = coef;
        if (n > m)
            p1 = p0 * x * (2*m+1);
            if (n > m+1)
                p = p1;
                pminus = p0;
                for i = 2:n1
                    pplus = 2.0*x*p - pminus + (2*m-1)*(x*p-pminus)/i;
                    pminus = p;
                    p = pplus;
                end
                pns = pplus;
            else
                pns = p1;
            end
        else
            pns = p0;
        end
    end
else
    pns = 0;
end


return
end

