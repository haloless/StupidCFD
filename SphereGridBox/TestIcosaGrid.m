
clear all;


[point,edge,face] = IcosahedralRegular();
edge = edge';
face = face';

npoint = size(point,2);
nedge = size(edge,1);
nface = size(face,1);

center = zeros(3,nface);
normvec = zeros(3,nface);
for i = 1:nface
    v1 = point(:,face(i,1));
    v2 = point(:,face(i,2));
    v3 = point(:,face(i,3));
    center(:,i) = (v1+v2+v3)/3;
    normvec(:,i) = cross(v2-v1,v3-v1); normvec(:,i) = normvec(:,i) ./ norm(normvec(:,i));
end

hfig = figure;
% trimesh(face, point(1,:)*0.9,point(2,:)*0.9,point(3,:)*0.9, ones(npoint,1));
axis equal;
axis([-1.1 1.1 -1.1 1.1 -1.1 1.1]);
hold on;
% quiver3(center(1,:),center(2,:),center(3,:), normvec(1,:),normvec(2,:),normvec(3,:));
hold off;

% face->edge
face2edge = zeros(nface,3);
edge2face = zeros(nedge,2);
for iface = 1:nface
    i1 = face(iface,1);
    i2 = face(iface,2);
    i3 = face(iface,3);
    
    for jedge = 1:nedge 
        ja = edge(jedge,1);
        jb = edge(jedge,2);
        
        % (i1,i2) 
        if i1==ja && i2==jb
            face2edge(iface,1) = jedge;
            edge2face(jedge,1) = iface;
        elseif i1==jb && i2==ja
            face2edge(iface,1) = -jedge;
            edge2face(jedge,2) = iface;
        end
        
        % i2,i3
        if i2==ja && i3==jb
            face2edge(iface,2) = jedge;
            edge2face(jedge,1) = iface;
        elseif i2==jb && i3==ja
            face2edge(iface,2) = -jedge;
            edge2face(jedge,2) = iface;
        end
        
        % i3,i1
        if i3==ja && i1==jb
            face2edge(iface,3) = jedge;
            edge2face(jedge,1) = iface;
        elseif i3==jb && i1==ja
            face2edge(iface,3) = -jedge;
            edge2face(jedge,2) = iface;
        end
    end
end

%
nref = 5;

%
nnode = npoint;
node = point;
edgenode = zeros(nedge,nref-1);
% refine edge
for iedge = 1:nedge
    ia = edge(iedge,1);
    ib = edge(iedge,2);
    xa = point(:,ia);
    xb = point(:,ib);
    
    for jref = 1:nref-1
        frac = jref / nref;
        xfrac = (1-frac)*xa + frac*xb;
        xfrac = xfrac ./ norm(xfrac);
        
        nnode = nnode + 1;
        node(:,nnode) = xfrac;
        edgenode(iedge,jref) = nnode;
    end
end

%
% facenode = zeros(nface,(nref+2)*(nref+1)/2);
elem = [];
nelem = 0;
%
for iface = 1:nface
    ia = face(iface,1);
    ib = face(iface,2);
    ic = face(iface,3);
    xa = point(:,ia);
    xb = point(:,ib);
    xc = point(:,ic);
    
    facenode = zeros(nref+1,nref+1);
    
    % cnt = 0;
    for v = 0:nref
    for u = 0:nref-v
        w = nref - u - v;
        % cnt = cnt + 1;
        
        
        if u==0 && v==0
            facenode(u+1,v+1) = ia;
        elseif v==0 && w==0
            facenode(u+1,v+1) = ib;
        elseif w==0 && u==0
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
        else % inside 
            nnode = nnode + 1;
            xnode = (w*xa + u*xb + v*xc)/nref;
            xnode = xnode / norm(xnode);
            
            facenode(u+1,v+1) = nnode;
            node(:,nnode) = xnode;
        end
    end
    end
    
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


if 0
    figure(hfig);
    hold on;
    plot3(node(1,:),node(2,:),node(3,:),'x');
    for v = 0:nref
    for u = 0:nref-v
        inode = facenode(u+1,v+1);
        text(node(1,inode),node(2,inode),node(3,inode),int2str(inode));
    end
    end
    hold off;
end

if 1
    figure(hfig);
    hold on;
    trimesh(elem,node(1,:),node(2,:),node(3,:),zeros(nnode,1));
    hold off;
end

center = zeros(3,nelem);
normvec = zeros(3,nelem);
for i = 1:nelem
    v1 = node(:,elem(i,1));
    v2 = node(:,elem(i,2));
    v3 = node(:,elem(i,3));
    center(:,i) = (v1+v2+v3)/3;
    normvec(:,i) = cross(v2-v1,v3-v1); normvec(:,i) = normvec(:,i) ./ norm(normvec(:,i));
end

if 1
    figure(hfig);
    hold on;
    quiver3(center(1,:),center(2,:),center(3,:), normvec(1,:),normvec(2,:),normvec(3,:));
    hold off;
end






