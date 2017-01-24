function [np] = IcosahedralPointsNum(factor)
%
% 12 vertex, 30 edge, 20 face
% 

np = 12 + 10*3*(factor-1) + 10*(factor-2)*(factor-1);

return
end

