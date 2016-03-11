
function [ alpha ] = StdChopGetIntercept (m1, m2, m3, vol)

CubeChoppingGlobals;

m12 = m1 + m2;
% mt = min(m12, m3);

v1 = m1^2 / max(6*m2*m3,ChopEps);
if (0<=vol && vol<v1)
    alpha = (6*m1*m2*m3*vol)^(1.0/3.0);
else
    v2 = v1 + (m2-m1)/(2*m3);
    if (v1<=vol && vol<v2)
        alpha = 0.5 * (m1 + sqrt(m1^2+8*m2*m3*(vol-v1)));
    else
        if (m3 < m12)
            v3 = (m3^2*(3*m12-m3)+m1^2*(m1-3*m3)+m2^2*(m2-3*m3)) / (6*m1*m2*m3);
        else
            v3 = m12 / (2*m3);
        end
        if (v2<=vol && vol<v3) 
            a3 = -1;
            a2 = 3*m12;
            a1 = -3*(m1^2+m2^2);
            a0 = m1^3 + m2^3 - 6*m1*m2*m3*vol;
            alpha = QubicSol(a3,a2,a1,a0);
        else
            if (v3<=vol && vol<=0.5)
                if (m3 < m12)
                    a3 = -2;
                    a2 = 3;
                    a1 = -3 * (m1^2+m2^2+m3^2);
                    a0 = m1^3+m2^3+m3^3 - 6*m1*m2*m3*vol;
                    alpha = QubicSol(a3,a2,a1,a0);
                else 
                    alpha = m3*vol + 0.5*m12;
                end
            else % vol>0.5
                alpha = StdChopGetIntercept(m1,m2,m3, 1.0-vol);
                alpha = 1.0-alpha;
            end
        end
    end
end

return
end


function [ sol ] = QubicSol(a3,a2,a1,a0)
a2 = a2 / a3;
a1 = a1 / a3;
a0 = a0 / a3;
p0 = a1/3 - (a2/3)^2;
q0 = (a1*a2-3*a0)/6 - (a2/3)^3;
theta = 1.0/3.0 * acos(q0/sqrt(-p0^3));
sol = sqrt(-p0) * (sqrt(3)*sin(theta)-cos(theta)) - a2/3;
end

