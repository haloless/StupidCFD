
function [u,v,w] = VortVel(x,y,z,t)

ProblemGlobals;
% U = 1.0;
% T = 8.0;

sx = sin(pi.*x);
sy = sin(pi.*y);
cx = cos(pi.*x);
cy = cos(pi.*y);

u = -U*cos(pi*t/T) .* (sx.^2) .* (2*sy.*cy);
v =  U*cos(pi*t/T) .* (sy.^2) .* (2*sx.*cx);
w = zeros(size(u));

return
end
