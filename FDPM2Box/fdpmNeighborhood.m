
function [ neigh,dist,rx,ry,re ] = fdpmNeighborhood (target,nodes,re)
%fdpmNeighborhood
% Get
% If $TARGET is an integer, then used as an index in $NODES
% otherwise used as coordinate

% check input target point
% if a scalar integer, then it's an index in the nodes array
% otherwise it's coordinate
if isscalar(target)
    isnode = true;
else
    isnode = false;
end

if isnode
    xc = nodes(target,1);
    yc = nodes(target,2);
else
    xc = target(1);
    yc = target(2);
end

rx = nodes(:,1)-xc;
ry = nodes(:,2)-yc;
dist = sqrt(rx.^2 + ry.^2);

if (~isscalar(re))
    % re is a vector, to be symmetrized
    rc = re(target);
    re = 1/2 * (re + rc);
end

flag = (dist < re);
if isnode
    % exclude target particle itself
    flag(target) = 0;
end

neigh = find(flag);

return
end
