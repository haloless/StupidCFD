
clear all;


hfig = figure;
% trimesh(face, point(1,:)*0.9,point(2,:)*0.9,point(3,:)*0.9, ones(npoint,1));
axis equal;
axis([-1.1 1.1 -1.1 1.1 -1.1 1.1]);
hold on;
% quiver3(center(1,:),center(2,:),center(3,:), normvec(1,:),normvec(2,:),normvec(3,:));
hold off;


%
% nref = 1;
nref = 2;
% nref = 3;

[node,elem] = IcosahedralGrid(nref);
nnode = size(node,2);
nelem = size(elem,1);

% kstiff = 1.0;
kstiff = 0.01;
vdamp = 0.01;
vel = zeros(size(node));
vel = 1.0e-2 * (rand(size(node))-1);
acc = zeros(size(node));
dt = 1.0e-3;
for step = 1:100000
    acc(:) = 0;
    for ielem = 1:nelem
        for ii = 1:3
            i = elem(ielem, ii);
            if ii == 3
                j = elem(ielem, 1);
            else
                j = elem(ielem, ii+1);
            end
            if (i < j)
                xij = node(:,i) - node(:,j);
                dij = norm(xij);
                nij = xij ./ dij;
                vij = vel(:,i) - vel(:,j);
                
                fij = kstiff / dij^2 * nij - vdamp*dot(vij,nij)*nij;
                acc(:,i) = acc(:,i) + fij;
                acc(:,j) = acc(:,j) - fij;
            end
        end
    end
    
    vel = vel + dt.*acc;
    node = node + dt.*vel;
    
    % project to sphere
    for i = 1:nnode
        ri = norm(node(:,i));
        node(:,i) = node(:,i) / ri;
    end
    
    if mod(step,100) == 0
        figure(hfig);
        % hold on;
        trimesh(elem,node(1,:),node(2,:),node(3,:),zeros(nnode,1));
        % hold off;
        title(['step=',int2str(step)]);
        drawnow;
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





