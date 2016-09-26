
function [ vol ] = StdChopGetVolume2D(m1,m2,alpha)
% Core function in 2D. Require m1+m2=1; m1<=m2

if alpha<0
    vol = 0.0;
elseif 0<=alpha & alpha<m1
    vol = alpha^2 / (2*m1*m2);
elseif m1<=alpha & alpha<=0.5
    vol = alpha/m2 - m1/(2*m2);
else % alpha>0.5
    % use complement
    vol = StdChopGetVolume2D(m1,m2,1-alpha);
    vol = 1.0 - vol;
return
end

