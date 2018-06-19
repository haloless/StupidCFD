function [node,elem] = IcosahedralGrid(nref)
% Function return sphere grid by refining icosahedron.
% nref=1, original 20-face regular icosahedron
%

% load regular icosahedron
[point,~,face] = IcosahedralRegular();
face = face';

% 
npoint = size(point,2);
nface = size(face,1);

% 
[edge,face2edge,edge2face] = MeshConnection(face);
nedge = size(edge,1);


% nodes 
nnode = npoint;
node = point;
% elements
nelem = 0;
elem = [];

map_to_sphere = 1;


%
% 1st pass, refine edges
% add new points and number them
%

edgenode = zeros(nedge,nref-1);

for iedge = 1:nedge
    ia = edge(iedge,1);
    ib = edge(iedge,2);
    xa = point(:,ia);
    xb = point(:,ib);
    
    for jref = 1:nref-1
        frac = jref / nref;
        xfrac = (1-frac)*xa + frac*xb;
        if map_to_sphere
            xfrac = xfrac ./ norm(xfrac);
        end
        
        nnode = nnode + 1;
        node(:,nnode) = xfrac;
        edgenode(iedge,jref) = nnode;
    end
end

%
% 2nd pass, refine faces
% add new points that are not on the face borders and number them
% then create elements
%
for iface = 1:nface
    ia = face(iface,1);
    ib = face(iface,2);
    ic = face(iface,3);
    xa = point(:,ia);
    xb = point(:,ib);
    xc = point(:,ic);
    
    facenode = zeros(nref+1,nref+1);
    
	% save node IDs on the barycentric coordinate of current face
    for v = 0:nref
    for u = 0:nref-v
        w = nref - u - v;        
        
        if u==0 && v==0 % point 1
            facenode(u+1,v+1) = ia;
        elseif v==0 && w==0 % point 2
            facenode(u+1,v+1) = ib;
        elseif w==0 && u==0 % point 3
            facenode(u+1,v+1) = ic;
        elseif v==0 % edge 12
            jedge = face2edge(iface,1);
            if jedge > 0
                facenode(u+1,v+1) = edgenode(jedge,u);
            else
                facenode(u+1,v+1) = edgenode(-jedge,nref-u);
            end
        elseif w==0 % edge 23
            jedge = face2edge(iface,2);
            if jedge > 0
                facenode(u+1,v+1) = edgenode(jedge,v);
            else
                facenode(u+1,v+1) = edgenode(-jedge,nref-v);
            end
        elseif u==0 % edge 31
            jedge = face2edge(iface,3);
            if jedge > 0
                facenode(u+1,v+1) = edgenode(jedge,w);
            else
                facenode(u+1,v+1) = edgenode(-jedge,nref-w);
            end
        else 
			% inside, the node has to be created
            nnode = nnode + 1;
            xnode = (w*xa + u*xb + v*xc)/nref;
            if map_to_sphere
                xnode = xnode / norm(xnode);
            end
            
            node(:,nnode) = xnode;
            facenode(u+1,v+1) = nnode;
        end
    end
    end
    
	% create sub-elements
    for v = 0:nref-1
    for u = 0:nref-1-v
        i00 = facenode(u+1,v+1);
        i10 = facenode(u+2,v+1);
        i01 = facenode(u+1,v+2);
        i11 = facenode(u+2,v+2);
        
        nelem = nelem + 1;
        elem(nelem,:) = [i00,i10,i01];
        
        if u+1+v+1 <= nref
            nelem = nelem + 1;
            elem(nelem,:) = [i10,i11,i01];
        end
    end
    end
end



return
end



