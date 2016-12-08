
function [n] = SpiralPointsEstimNum(dist)
% Estimate for unit sphere (a=1)
% The typical distance is given relative to a.

coef = 2.0 / sqrt(3.0);
n = ceil(coef * 4*pi / dist^2);


return
end

