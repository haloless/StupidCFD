function [sdf] = MakeGridSDF(shape,node,elem)


lx = shape.a * 6;
ly = shape.b * 6;
xmin = shape.xc - lx/2;
xmax = shape.xc + lx/2;
ymin = shape.yc - ly/2;
ymax = shape.yc + ly/2;

dh = min(shape.a,shape.b) / 20;
[xg,yg] = ndgrid(xmin:dh:xmax,ymin:dh:ymax);
nx = size(xg,1);
ny = size(xg,2);

np = size(node,2);
ne = size(elem,1);

centelem = zeros(2,ne);
for i = 1:ne
    i1 = elem(i,1);
    i2 = elem(i,2);
    centelem(:,i) = (node(:,i1)+node(:,i2))/2;
end

nvecelem = zeros(2,ne);
for i = 1:ne
    i1 = elem(i,1);
    i2 = elem(i,2);
    xx = node(1,i2) - node(1,i1);
    yy = node(2,i2) - node(2,i1);
    nn = [yy;-xx];
    nvecelem(:,i) = nn ./ norm(nn);
end

nvecnode = zeros(2,np);
for i = 1:ne
    i1 = elem(i,1);
    i2 = elem(i,2);
    nvecnode(:,i1) = nvecnode(:,i1) + nvecelem(:,i);
    nvecnode(:,i2) = nvecnode(:,i2) + nvecelem(:,i);
end
for i = 1:np
    nvecnode(:,i1) = nvecnode(:,i1) ./ norm(nvecnode(:,i1));
end

if 0
    hold on;
    quiver(node(1,:),node(2,:),nvecnode(1,:),nvecnode(2,:));
    quiver(centelem(1,:),centelem(2,:), nvecelem(1,:),nvecelem(2,:));
    hold off;
end

phig = zeros(nx,ny);
far = 99999;
phig(:,:) = far;

% element distance
for i = 1:nx
for j = 1:ny
    pos = [xg(i,j);yg(i,j)];
    
    for kelem = 1:ne
        vcen = centelem(:,kelem);
        ncen = nvecelem(:,kelem);
        
        k1 = elem(kelem,1);
        k2 = elem(kelem,2);
        v1 = node(:,k1);
        v2 = node(:,k2);
        
        
        vrel = pos - vcen;
        vproj = dot(vrel,ncen) * ncen;
        vpara = vrel - vproj;
        
        if norm(vpara) > 0.5*norm(v2-v1)
            dist1 = norm(pos-v1);
            dist2 = norm(pos-v2);
            if dist1 < dist2
                dist = sign(dot(pos-v1,nvecnode(:,k1))) * dist1;
            else
                dist = sign(dot(pos-v2,nvecnode(:,k2))) * dist2;
            end
        else
            dist = norm(vproj) * sign(dot(vproj,ncen));
        end
        
        if abs(phig(i,j)) > abs(dist)
            phig(i,j) = dist;
        end
    end
end
end



% phig = zeros(nx,ny);
% for i = 1:nx
% for j = 1:ny
	% phig(i,j) = ShapePotential(shape,xg(i,j),yg(i,j));
% end
% end

if 1
    % reinitialize to get SDF
    band = lx * 0.2;
    phig = ImplicitFuncReinit(phig,nx,ny,dh,dh, band);
end

sdf = struct();
sdf.xmin = xmin;
sdf.xmax = xmax;
sdf.ymin = ymin;
sdf.ymax = ymax;
sdf.dh = dh;
sdf.nx = nx;
sdf.ny = ny;
sdf.xg = xg;
sdf.yg = yg;
sdf.phig = phig;



return
end

