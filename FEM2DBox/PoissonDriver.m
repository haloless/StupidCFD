% ## Copyright (C) 2013 homu
% ## 
% ## This program is free software; you can redistribute it and/or modify
% ## it under the terms of the GNU General Public License as published by
% ## the Free Software Foundation; either version 3 of the License, or
% ## (at your option) any later version.
% ## 
% ## This program is distributed in the hope that it will be useful,
% ## but WITHOUT ANY WARRANTY; without even the implied warranty of
% ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% ## GNU General Public License for more details.
% ## 
% ## You should have received a copy of the GNU General Public License
% ## along with Octave; see the file COPYING.  If not, see
% ## <http://www.gnu.org/licenses/>.

% ## BeamDriver

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-30

% clc;
clear;


%% source function

srctype = 'sinsin';
% srctype = 'linear';
switch srctype
case 'sinsin' 
    % f=sin(m*pi*x) * sin(m*pi*y)
    mm = 1; nn = 1;
    % mm = 2; nn = 3;
    fun_fext = @(x,y) (mm^2+nn^2)/2 * sin(mm*pi*x) .* sin(nn*pi*y);
    fun_phi = @(x,y) 1/(2*pi^2) * sin(mm*pi*x) .* sin(nn*pi*y);
    fun_phix = @(x,y) mm/(2*pi) * cos(mm*pi*x) .* sin(nn*pi*y);
    fun_phiy = @(x,y) nn/(2*pi) * sin(mm*pi*x) .* cos(nn*pi*y);
case 'linear'
    % f = a*x + b*y + c
    aa = 1; bb = 3; cc = 0;
    fun_fext = @(x,y) zeros(size(x));
    fun_phi = @(x,y) cc + aa*x + bb*y;
    fun_phix = @(x,y) aa * ones(size(x));
    fun_phiy = @(x,y) bb * ones(size(y));
otherwise
    error('unknown srctype: %s', srctype);
end



%% 
Lx = 1;
Ly = 1;

% nx = 10;
% ny = 10;
% nx = 20;
% ny = 20;
nx = 30;
ny = 30;

nx1 = nx + 1;
ny1 = ny + 1;


%% 
elem_type = 'Q4';
elem_numNodes = 4;
elem_numGauss = 4;


%% generate mesh
numElems = nx * ny;
numNodes = (nx+1) * (ny+1);
[meshX,meshY] = ndgrid(linspace(0,Lx,nx+1),linspace(0,Ly,ny+1));
nodeX = reshape(meshX,[],1);
nodeY = reshape(meshY,[],1);
nodes = [nodeX, nodeY];
% DOF with [phi]
numDofs = numNodes;

% build connectivity
T = zeros(numElems,elem_numNodes);
iElem = 0;
for j = 1:ny
    for i = 1:nx
        iElem = iElem + 1;
        n1 = i + (j-1)*(nx+1);
        T(iElem,:) = [n1, n1+1, n1+1+(nx+1), n1+nx+1];
    end
end

if (1)
    figure; plot_mesh(nodes,T,elem_type);
    % return
end

% analytical solution 

phis = fun_phi(nodeX,nodeY);


%% global stiffness matrix
Kmat = sparse(numDofs,numDofs);
rhs = zeros(numDofs,1);

disp(['Build global stiffness matrix K']);
for iElem = 1:numElems
    iNodeIds = T(iElem,:);
    iNodePos = nodes(iNodeIds,:);
    iDofs = iNodeIds';
    
    Ke = zeros(elem_numNodes,elem_numNodes);
    fe = zeros(elem_numNodes,1);
    
    % load Gauss points
    [gpoints,weights] = GaussPoint(elem_numNodes,elem_numGauss);
    % axi = gpoints(:,1);
    % aeta = gpoints(:,2);
    
    for kGauss = 1:elem_numGauss
        localCoord = gpoints(kGauss,:);
        xi = gpoints(kGauss,1);
        eta = gpoints(kGauss,2);
        wgt = weights(kGauss);
        
        % shape function in local world
        [N,dN] = LagrangeBasis(elem_type,localCoord);
        
        % Jacobian matrix
        [Jmat,invJ,detJ] = JacobianMatrix(dN,iNodePos);
        
        % shape function in real world
        dN = (Jmat \ dN')';
        dNdX = dN(:,1);
        dNdY = dN(:,2);
        
        %
        Ke = Ke + (dN*dN') * detJ*wgt;
        
        % external source
        xg = dot(N,iNodePos(:,1));
        yg = dot(N,iNodePos(:,2));
        fg = fun_fext(xg,yg);
        
        fe = fe + (N) * fg * detJ*wgt;
    end
    
    % assemble
    Kmat(iDofs,iDofs) = Kmat(iDofs,iDofs) + Ke;
    
    rhs(iDofs) = rhs(iDofs) + fe;
end

disp(['Apply BC.']);

tol = 1.0e-6;
disp_ids = find(abs(nodeX-0)<tol | abs(nodeX-Lx)<tol | abs(nodeY-0)<tol | abs(nodeY-Ly)<tol);
disp_dofs = disp_ids(:);
disp_vals = phis(disp_dofs);

trac_ids = [];
trac_dofs = [];
trac_vals = [];
if (1)
    figure;
    plot(nodeX(disp_ids),nodeY(disp_ids),'x', nodeX(trac_ids),nodeY(trac_ids),'o'); 
    legend('disp','trac');
    axis equal;
end

free_dofs = ones(numDofs,1)';
free_dofs(disp_dofs) = 0;
free_dofs = find(free_dofs);


% set natural BC
rhs(trac_dofs) = trac_vals;

% enforce essential BC
rhs(free_dofs) = rhs(free_dofs) - Kmat(free_dofs,disp_dofs)*disp_vals;

sol = zeros(numDofs,1);

sol(free_dofs) = Kmat(free_dofs,free_dofs) \ rhs(free_dofs);
sol(disp_dofs) = disp_vals;


disp(['Solving equations.']);
% sol = Kmat \ rhs;
% sol = reshape(sol,2,[]);
% fem_u = sol(1,:)';
% fem_v = sol(2,:)';

% fint = Kmat * sol(:);
% fint = reshape(fint, 2,[]);
% fintx = reshape(fint(1,:),nx+1,ny+1);
% finty = reshape(fint(2,:),nx+1,ny+1);


if (1)
    xx = reshape(nodeX, nx1,ny1);
    yy = reshape(nodeY, nx1,ny1);
    
    figure;
    zz = reshape(sol, nx1,ny1);
    surf(xx,yy,zz,'FaceColor','interp');
    
    % plot(nodeX+fem_u,nodeY+fem_v,'o', nodeX+disp_u,nodeY+disp_v,'x'); 
    % legend('FEM','patch'); axis equal;
    
    figure;
    surf(xx,yy,reshape(sol-phis, nx1,ny1));
    % subplot(2,1,1);
    % surf(meshX',meshY',reshape(fem_u-disp_u,nx+1,ny+1)');
    % title('err-u');
    % subplot(2,1,2);
    % surf(meshX',meshY',reshape(fem_v-disp_v,nx+1,ny+1)');
    % title('err-v');
    % axis equal;
end


if 1
    % check stiffness K matrix eigen structure
    
    figure;
    
    % we copy the K matrix and set Dirichlet BC for boundaries
    Ktan = Kmat;
    Ktan(sub2ind([numDofs,numDofs],disp_dofs,disp_dofs)) = 1.0e15;
    
    neig = 20;
    % [v,d] = eigs(Ktan, neig);
    [v,d] = eigs(Ktan, neig, 'sm'); % from min eigenvalue
    for ieig = 1:neig
        vv = v(:,ieig);
        surf(xx,yy, reshape(vv,nx1,ny1));
        title(['eig(',int2str(ieig),')=',num2str(d(ieig,ieig))]);
        pause;
    end
end



