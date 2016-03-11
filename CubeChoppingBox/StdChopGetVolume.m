
function [ vol ] = StdChopGetVolume (m1, m2, m3, alpha)

CubeChoppingGlobals;

m12 = m1 + m2;
mt = min(m12, m3);

if (alpha < 0)
    vol = 0;
elseif (0<=alpha && alpha<m1)
    vol = alpha^3 / (6*m1*m2*m3);
elseif (m1<=alpha && alpha<m2)
    v1 = m1^2 / max(6*m2*m3,ChopEps);
    vol = alpha*(alpha-m1)/(2*m2*m3) + v1;
elseif (m2<=alpha && alpha<mt)
    vol = alpha^2*(3*m12-alpha) + m1^2*(m1-3*alpha) + m2^2*(m2-3*alpha);
    vol = vol / (6*m1*m2*m3);
elseif (mt<=alpha && alpha<=0.5)
    if (m3 < m12) % mt=min(m3,m12)=m3
        vol = alpha^2*(3-2*alpha);
        vol = vol + m1^2*(m1-3*alpha) + m2^2*(m2-3*alpha) + m3^2*(m3-3*alpha);
        vol = vol / (6*m1*m2*m3);
    else % mt=min(m3,m12)=m12
        vol = (2*alpha-m12) / (2*m3);
    end
else % alpha>0.5
    % use its complement 
    vol = StdChopGetVolume(m1,m2,m3, 1.0-alpha);
    vol = 1.0 - vol;
end

return
end

