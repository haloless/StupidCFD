function [tri] = SpherePointsToGrid(x,y,z)

% this is an easy hack
% if input points are all on sphere surface,
% then their convex hull is equivalent to a sphere Delaunay
tri = convhull(x,y,z);


return
end

