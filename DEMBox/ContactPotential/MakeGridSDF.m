function [sdf] = MakeGridSDF(shape, node,elem, do_reinit)


lx = shape.a * 6;
ly = shape.b * 6;
xmin = shape.xc - lx/2;
xmax = shape.xc + lx/2;
ymin = shape.yc - ly/2;
ymax = shape.yc + ly/2;

% len = max(shape.a,shape.b);
dh = min(shape.a,shape.b) / 20;
[xg,yg] = ndgrid(xmin:dh:xmax,ymin:dh:ymax);
nx = size(xg,1);
ny = size(xg,2);

np = size(node,2);
ne = size(elem,1);

centelem = zeros(2,ne);
lenelem = zeros(ne,1);
for i = 1:ne
    i1 = elem(i,1);
    i2 = elem(i,2);
    centelem(:,i) = (node(:,i1)+node(:,i2))/2;
	lenelem(i) = norm(node(:,i1)-node(:,i2));
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

if 0
	% full distance
	for i = 1:nx
	for j = 1:ny
		pos = [xg(i,j);yg(i,j)];
		
		for kelem = 1:ne
			vcen = centelem(:,kelem);
			ncen = nvecelem(:,kelem);
			lcen = lenelem(kelem);
			
			k1 = elem(kelem,1);
			k2 = elem(kelem,2);
			v1 = node(:,k1);
			v2 = node(:,k2);
			
			
			vrel = pos - vcen;
			vproj = dot(vrel,ncen) * ncen;
			vpara = vrel - vproj;
			
			if norm(vpara) > 0.5*lcen
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
end

if 1
	% near distance
	len = 0.5 * max(shape.a,shape.b);
	margin = floor(len/dh);
	disp(['margin=',int2str(margin)]);
	
	for kelem = 1:ne
		
		vcen = centelem(:,kelem);
		ncen = nvecelem(:,kelem);
		lcen = lenelem(kelem);
		
		k1 = elem(kelem,1);
		k2 = elem(kelem,2);
		v1 = node(:,k1);
		v2 = node(:,k2);
		
		imin = floor((min(v1(1),v2(1))-xmin)/dh) + 1 - margin;
		jmin = floor((min(v1(2),v2(2))-ymin)/dh) + 1 - margin;
		imax = ceil((max(v1(1),v2(1))-xmin)/dh) + 1 + margin;
		jmax = ceil((max(v1(2),v2(2))-ymin)/dh) + 1 + margin;
		imin = max(imin, 1);
		jmin = max(jmin, 1);
		imax = min(imax, nx);
		jmax = min(jmax, ny);
		
		for i = imin:imax
		for j = jmin:jmax
			pos = [ xg(i,j); yg(i,j) ];
			
			vrel = pos - vcen;
			vproj = dot(vrel,ncen) * ncen;
			vpara = vrel - vproj;
			
			if norm(vpara) > 0.5*lcen
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
		
		if mod(kelem,10) == 0
			fprintf('elem=%d in %d\n',kelem,ne);
		end
	end
	
	% fill far distance
	
	tag = ones(nx,ny);
	tag(phig==far) = 0;
	
	chk = find(tag==0);
	
	iter = 0;
	while ~isempty(chk)
		for kchk = 1:numel(chk)
			ind = chk(kchk);
			[ic,jc] = ind2sub([nx,ny],ind);
			
			ok = 0;
			if ic > 1
				if tag(ic-1,jc) ~= 0
					phig(ic,jc) = phig(ic-1,jc);
					ok = 1;
				end
			end
			if ic < nx
				if tag(ic+1,jc) ~= 0
					phig(ic,jc) = phig(ic+1,jc);
					ok = 1;
				end
			end
			if jc > 1
				if tag(ic,jc-1) ~= 0
					phig(ic,jc) = phig(ic,jc-1);
					ok = 1;
				end
			end
			if jc < ny
				if tag(ic,jc+1) ~= 0
					phig(ic,jc) = phig(ic,jc+1);
					ok = 1;
				end
			end
			
			if ok
				tag(ic,jc) = 1;
				chk(kchk) = 0;
			end
		end
		
		chk = chk(chk>0);
		
		iter = iter + 1;
		if mod(iter,20)==0
			disp(['iter=',int2str(iter)]);
		end
	end
	
	
end


% reinitialize to get regularized SDF
if ~exist('do_reinit','var')
	do_reinit = 1;
end
if do_reinit
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

