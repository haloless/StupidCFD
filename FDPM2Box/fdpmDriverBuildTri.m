%%fdpmDriverBuildTri: Build triangulation of input nodes.
% A (numEdges, edgePair) can be provided to specify the boundary.
% Results
% * deltri: the Delaunay Triangulation object
% * numTri, tri: the triangle table
% * trep: the triangulation representation (consider boundary)
% * nodeArea: area of each node by dividing triangles
%

show_tri_plot = 1;


% perform triangulation
if ~exist('edgePair','var')
    warning('edgePair is empty, triangulation may not be valid');
    edgePair = [];
    numEdges = 0;
end

if isempty(edgePair)
    % triangulation without edge constraints
    deltri = DelaunayTri(nodePos);
    
    tri = deltri(:,:);
else
    % triangulation with edge constraints
    deltri = DelaunayTri(nodePos, edgePair);
    
    % judge in/out
    io = inOutStatus(deltri);
    
    % extract valid triangles 
    tri = deltri(io,:);
end

% number of triangles
numTri = size(tri,1);

% the triangulation representation object
trep = TriRep(tri, nodeX,nodeY);

% calculate 
% cc = circumcenters(trep);




if show_tri_plot
    figure;
    plot(nodeX,nodeY,'.');
    axis equal;
    
    hold on;
    
    if ~isempty(edgePair)
        for iedge = 1:numEdges
            plot(nodeX(edgePair(iedge,:)),nodeY(edgePair(iedge,:)), 'r','LineWidth',2.5);
        end
    end
    
    triplot(tri, nodeX,nodeY);
    
    % plot(cc(:,1),cc(:,2),'g.');
    
    hold off;
    
    % pause;
end


if 1
    % calculate area associated with each node
    nodeArea = zeros(numNodes,1);
    
    for i = 1:numTri
        i1 = tri(i,1); 
        i2 = tri(i,2); 
        i3 = tri(i,3);
        
        l1 = norm(nodePos(i2,:)-nodePos(i3,:));
        l2 = norm(nodePos(i3,:)-nodePos(i1,:));
        l3 = norm(nodePos(i1,:)-nodePos(i2,:));
        
        ll = [l1,l2,l3];
        
        ii = 1;
        if l2 > ll(ii)
            ii = 2;
        end
        if l3 > ll(ii)
            ii = 3;
        end
        
        if ii == 1
            ia = i1; ib = i2; ic = i3;
            la = l1; lb = l2; lc = l3;
        elseif ii == 2
            ia = i2; ib = i3; ic = i1;
            la = l2; lb = l3; lc = l1;
        else
            ia = i3; ib = i1; ic = i2;
            la = l3; lb = l1; lc = l2;
        end
        assert(la>=lb && la>=lc);
        
        ta = acos((lb^2+lc^2-la^2)/(2*lb*lc));
        if ta >= pi/2
            ss = 0.5*sin(ta)*lb*lc;
            nodeArea(ia) = nodeArea(ia) + ss/2;
            nodeArea(ib) = nodeArea(ib) + ss/4;
            nodeArea(ic) = nodeArea(ic) + ss/4;
        else
            tb = acos((lc^2+la^2-lb^2)/(2*lc*la));
            tc = pi - ta - tb;
            
            sa = la/2/tan(ta) * la/2;
            sb = lb/2/tan(tb) * lb/2;
            sc = lc/2/tan(tc) * lc/2;
            
            nodeArea(ia) = nodeArea(ia) + (sb+sc)/2;
            nodeArea(ib) = nodeArea(ib) + (sc+sa)/2;
            nodeArea(ic) = nodeArea(ic) + (sa+sb)/2;
        end
    end
end

if 1
    % in 2D, set volume = area
    % this may be corrected in following operations
    nodeVol = nodeArea;
end






