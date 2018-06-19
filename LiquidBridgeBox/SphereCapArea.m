function [area] = SphereCapArea(a,h)
% Area of spherical cap.
% a: radius of cap (not sphere radius)
% h: height of cap

r = (a^2+h^2)/(h*2);
area = 2*pi * r*h;

return
end


