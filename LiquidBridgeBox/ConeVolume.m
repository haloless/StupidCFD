function [V] = ConeVolume(r1,r2,h)
% Volume of Frustum Cone with (r1,r2,h)

V = pi/3 .* h .* (r1.^2 + r2.^2 + r1.*r2);
return
end

