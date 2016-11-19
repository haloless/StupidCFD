
clear all;

BemMeshGlobals;

% kappa = 10.0;
kappa = 1.0;
% kappa = 0.5;
% kappa = 0.2;
% kappa = 0.1;
% kappa = 0.05;
% kappa = 0.01;
% kappa = 0.001;

cutoff = 3.0;
% cutoff = 9999;

input_file = 'tmp.dat'; 
% input_file = 'input.dat'; 
% input_file = 'input_ring_M360.dat'; 
% input_file = 'input_square_M400q.dat'; 
% input_file = 'input_square_M400.dat'; 

%
BemLoadMesh;

if (0)
for i = 1:nelem
    if bc(1,i) == bc_dir
        bc(2,i) = 1.0;
    elseif bc(1,i) == bc_neu
        bc(2,i) = 0.0;
    end
end
end

% plot elements
if (1)
    BemPlotMesh;
end

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

if (1)
    figure;
    hold on;
    plot3(y(1,:),y(2,:),zeros(npoint,1),'.');
    plot3(x(1,:),x(2,:),uval,'x-');
    plot3(xfield(1,:),xfield(2,:),f,'o');
    hold off;
    % axis equal;
end

if (0)
    % analytical solution 1D planar
    xlo = 0.0;
    xhi = 1.0;
    ulo = 1.0;
    uhi = 1.0;
    coef = [exp(kappa*xlo),exp(-kappa*xlo);exp(kappa*xhi),exp(-kappa*xhi)] \ [ulo;uhi];
    acoef = coef(1);
    bcoef = coef(2);
    xx = xfield(1,:)';
    fana = acoef.*exp(kappa.*xx) + bcoef.*exp(-kappa.*xx);
end
if (1)
    % analytical solution 1D radial
    rlo = 40;
    rhi = 80;
    ulo = 1.5;
    uhi = 1.0;
    coef = [besseli(0,kappa*rlo),besselk(0,kappa*rlo);besseli(0,kappa*rhi),besselk(0,kappa*rhi)];
    coef = coef \ [ulo;uhi];
    acoef = coef(1);
    bcoef = coef(2);
    xx = xfield(1,:)';
    fana = acoef.*besseli(0,kappa.*xx) + bcoef.*besselk(0,kappa.*xx);
end








