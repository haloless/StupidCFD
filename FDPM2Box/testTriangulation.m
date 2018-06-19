%testTriangulation
% 

clear;

Ra = 1;
Rb = 2;

nrad = 4;
nphi = 9;
% nrad = 8;
% nphi = 18;

drad = (Rb-Ra) / nrad;
dphi = pi/2 / nphi;

numNodes = 0;
nodeX = [];
nodeY = [];
% nodeInd = [];

nodeInd = @(irad,iphi) (irad+1)+(nrad+1)*iphi;

if 0
    % 1/4 cylinder, constraint on x- and y-axis
    for iphi = 0:nphi
    for irad = 0:nrad
        rr = Ra + irad * drad;
        tt = iphi * dphi;
        
        xx = rr * cos(tt);
        yy = rr * sin(tt);
        
        numNodes = numNodes + 1;
        nodeX(numNodes,1) = xx;
        nodeY(numNodes,1) = yy;
        % nodeInd(irad,iphi) = numNodes;
    end
    end
else
    % 1/4 cylinder, constraint on x- and y-axis
    for iphi = 0:nphi
    for irad = 0:nrad
        xx = irad * drad;
        yy = iphi * drad;
        
        numNodes = numNodes + 1;
        nodeX(numNodes,1) = xx;
        nodeY(numNodes,1) = yy;
        % nodeInd(irad,iphi) = numNodes;
    end
    end
end

% define edge
numEdges = 0;
edgePair = [];
for iphi = [0,nphi]
    for irad = 0:nrad-1
        numEdges = numEdges + 1;
        ind1 = nodeInd(irad,iphi);
        ind2 = nodeInd(irad+1,iphi);
        edgePair(numEdges,:) = [ind1,ind2];
    end
end
for irad = [0,nrad]
    for iphi = 0:nphi-1
        numEdges = numEdges + 1;
        ind1 = nodeInd(irad,iphi);
        ind2 = nodeInd(irad,iphi+1);
        edgePair(numEdges,:) = [ind1,ind2];
    end
end

nodePos = [nodeX,nodeY];

%
dt = DelaunayTri(nodePos,edgePair);

% judge in/out
io = inOutStatus(dt);

% 
tri = dt(io,:);

% 
trep = TriRep(tri, nodeX,nodeY);

cc = circumcenters(trep);

%
% dt.Triangulation = tri;

% [vo,ro] = voronoiDiagram(dt);


if 1
    figure;
    plot(nodeX,nodeY,'.');
    axis equal;
    
    hold on;
    
    for iedge = 1:numEdges
        plot(nodeX(edgePair(iedge,:)),nodeY(edgePair(iedge,:)), 'r','LineWidth',2.5);
    end
    
    % triplot(tri,nodeX,nodeY);
    % triplot(dt);
    triplot(tri, nodeX,nodeY);
    
    plot(cc(:,1),cc(:,2),'g.');
    
    hold off;
    
    % pause;
end

nodeArea = zeros(numNodes,1);

numTri = size(tri,1);
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








