
clear all;

BemMeshGlobals;

%
% Single sphere problem
% phi(r) = phi0*a*exp(ka) * exp(-kr)/r
% dphi/dr = -phi0*a*exp(ka) * exp(-kr) * (1/r^2 + k/r)
%


% kappa = 10.0;
% kappa = 2.0;
% kappa = 1.0;
kappa = 0.5;
% kappa = 0.2;
% kappa = 0.1;
% kappa = 0.05;
% kappa = 0.01;
% kappa = 0.001;

disp('Generate mesh');

% generate sphere surface
addpath('../../SphereGridBox');
if 1
	% use refined icosahedron
	vert0 = IcosahedralPoints(2);
	% vert0 = IcosahedralPoints(3);
	% vert0 = IcosahedralPoints(4);
	face0 = convhull(vert0(1,:),vert0(2,:),vert0(3,:));
end
% build surface grid
if 0
	vert = vert0;
	face = face0;
end
if 1
	v1 = [vert0(1,:)-1.5;vert0(2,:);vert0(3,:)];
	v2 = [vert0(1,:)+1.5;vert0(2,:);vert0(3,:)];
	f1 = face0;
	f2 = f1 + size(v1,2);
	vert = [v1, v2];
	face = [f1; f2];

end

% vertex/face number
nvert = size(vert,2);
nface = size(face,1);
disp(['vertex=',int2str(nvert)]);
disp(['face=',int2str(nface)]);

% BC
bc = zeros(2,nface);
bc(1,:) = 1;
bc(2,:) = 1.0;
% bc(1,:) = 2;
% bc(2,:) = -3;

% mesh status
facecent = zeros(3,nface);
facenvec = zeros(3,nface);
facearea = zeros(nface,1);
facelmax = zeros(nface,1);
facelmin = zeros(nface,1);
for i = 1:nface
	i1 = face(i,1);
	i2 = face(i,2);
	i3 = face(i,3);
	
	v1 = vert(:,i1);
	v2 = vert(:,i2);
	v3 = vert(:,i3);
	
	% centroid
	facecent(:,i) = (v1+v2+v3) / 3.0;
	
	% side length
	dl12 = norm(v2-v1);
	dl23 = norm(v3-v2);
	dl31 = norm(v1-v3);
	facelmax(i) = max([dl12,dl23,dl31]);
	facelmin(i) = min([dl12,dl23,dl31]);
	
	% area and normal
	ea = vert(:,i2) - vert(:,i1);
	eb = vert(:,i3) - vert(:,i1);
	ec = cross(ea,eb);
	facearea(i) = norm(ec) * 0.5;
	facenvec(:,i) = ec ./ norm(ec);
	
end



% build matrix
disp('Build matrix');
tic;
BemBuildConstElem;
% BemBuildNodePatch;
toc;

% solve
% sol = amat \ bvec;
sol = gmres(amat,bvec,20,1.0e-8,2000);

if 1
	figure;
	trimesh(face,vert(1,:),vert(2,:),vert(3,:));
	% trimesh(face,vert(1,:),vert(2,:),vert(3,:),sol);
	hold on;
	quiver3(facecent(1,:),facecent(2,:),facecent(3,:),facenvec(1,:),facenvec(2,:),facenvec(3,:));
	hold off;
	axis equal;
	xlabel('x');ylabel('y');zlabel('z');
end

return

%
% evaluate fields
%
disp('Fields');
xgrid = linspace(-2,2,21);
ygrid = linspace(-2,2,21);
[xgrid,ygrid] = ndgrid(xgrid,ygrid);
fgrid = zeros(size(xgrid));

tic;
for kelem = 1:nface
	if bc(1,kelem) == 1
		f0 = bc(2,kelem);
		df0 = sol(kelem);
	elseif bc(1,kelem) == 2
		f0 = sol(kelem);
		df0 = bc(2,kelem);
	end
	
	for i = 1:numel(fgrid)
		level = 1;
		[aa,bb] = BemCoef2([xgrid(i),ygrid(i),0], kelem, level);
		fgrid(i) = fgrid(i) + aa*df0 - bb*f0;
	end
	
	if mod(kelem,20) == 0
		disp(['elem=',int2str(kelem)]);
	end
end
toc;

figure;
contourf(xgrid,ygrid,fgrid);
axis equal;

return


%
% matrix and vector
%
disp('Build linear equations');
tic;
Amat = zeros(nelem,nelem);
bvec = zeros(nelem,1);
for j = 1:nelem
    for i = 1:nelem
        % 
        rij = norm(x(:,j)-x(:,i));
        is_singular = (rij <= lenmax*cutoff);
        
        [aa,bb] = BemCoef(x(1,i),x(2,i),j, is_singular);
        if i == j
            bb = 0.5;
        end
        
        if bc(1,j) == bc_dir
            Amat(i,j) = Amat(i,j) - aa;
            bvec(i) = bvec(i) - bb*bc(2,j);
        elseif bc(1,j) == bc_neu
            Amat(i,j) = Amat(i,j) + bb;
            bvec(i) = bvec(i) + aa*bc(2,j);
        end
    end
end
toc;

%
% solution
%
disp('Solver');
tic;
sol = Amat \ bvec;
u = sol;
toc;

% u contains both value and derivative
% replace all Dirichlet BC
uval = sol;
for i = 1:nelem
    if bc(1,i) == bc_dir
        uval(i) = bc(2,i);
    end
end



%
% evaluate fields
%
disp('Fields');
tic;
f = zeros(nfield,1);
for j = 1:nelem
    if bc(1,j) == bc_dir
        f0 = bc(2,j);
        df0 = u(j);
    elseif bc(1,j) == bc_neu
        f0 = u(j);
        df0 = bc(2,j);
    end

    for i = 1:nfield
        % 
        rij = norm(x(:,j)-xfield(:,i));
        is_singular = (rij <= lenmax*cutoff);
        
        [aa,bb] = BemCoef(xfield(1,i),xfield(2,i),j, is_singular);
        f(i) = f(i) + aa*df0 - bb*f0;
    end
end
toc;









