function [vol] = SphereCapVolume(a,h)
% Volume of spherical cap.
% a: radius of cap (not sphere radius)
% h: height of cap

vol = pi/6 * h * (3*a^2+h^2);

return
end


