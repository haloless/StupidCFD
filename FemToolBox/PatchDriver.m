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

%% 
E0 = 500;
nu0 = 0.25;

% plane stress
% Dmat = E0 / (1-nu0^2) * [1,nu0,0; nu0,1,0; 0,0,(1-nu0)/2];

% plane strain
Dmat = E0/(1+nu0)/(1-2*nu0) * [...
1-nu0, nu0,   0; 
nu0,   1-nu0, 0;
0,    0,    0.5-nu0];

%% 
Lx = 10;
Ly = 20;

nx = 10;
ny = 20;

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
% DOF with [u,v]
numDofs = numNodes*2;

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
end

% analytical solution of Timoshenko

% fun_u = @(x,y) 0.01*x;
% fun_v = @(x,y) 0.01*y;

fun_u = @(x,y) 0.1 + 0.1*x + 0.2*y;
fun_v = @(x,y) 0.05 + 0.15*x + 0.1*y;

disp_u = fun_u(nodeX,nodeY);
disp_v = fun_v(nodeX,nodeY);


%% global stiffness matrix
Kmat = sparse(numDofs,numDofs);
rhs = zeros(numDofs,1);

disp(['Build global stiffness matrix K']);
for iElem = 1:numElems
    iNodeIds = T(iElem,:);
    iNodePos = nodes(iNodeIds,:);
    iDofs = zipcol([iNodeIds*2-1;iNodeIds*2]');
    
    Ke = zeros(2*elem_numNodes,2*elem_numNodes);
    
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
        
        % B matrix
        Bmat = [zipcol(dNdX*[1,0])'; zipcol(dNdY*[0,1])'; zipcol([dNdY,dNdX])'];
        if (size(Bmat,1)~=3 | size(Bmat,2)~=2*elem_numNodes)
            error('Bmat error');
        end
        
        % BtDB
        Ke = Ke + (Bmat'*Dmat*Bmat) * detJ * wgt;
    end
    
    % assemble
    Kmat(iDofs,iDofs) = Kmat(iDofs,iDofs) + Ke;
end

disp(['Apply BC.']);

tol = 1.0e-6;
disp_ids = find(abs(nodeX-0)<tol | abs(nodeX-Lx)<tol | abs(nodeY-0)<tol | abs(nodeY-Ly)<tol);
disp_dofs = zipcol([disp_ids*2-1, disp_ids*2]);
disp_vals = zipcol([disp_u(disp_ids),disp_v(disp_ids)]);

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

% spI = speye(numNodes*2,numNodes*2);
% Kmat(disp_dofs,:) = spI(disp_dofs,:);
% rhs(disp_dofs) = disp_vals;

disp(['Solving equations.']);
% sol = Kmat \ rhs;
sol = reshape(sol,2,[]);
fem_u = sol(1,:)';
fem_v = sol(2,:)';

fint = Kmat * sol(:);
fint = reshape(fint, 2,[]);
fintx = reshape(fint(1,:),nx+1,ny+1);
finty = reshape(fint(2,:),nx+1,ny+1);


if (1)
    figure;
    plot(nodeX+fem_u,nodeY+fem_v,'o', nodeX+disp_u,nodeY+disp_v,'x'); 
    legend('FEM','patch'); axis equal;
    
    figure;
    subplot(2,1,1);
    surf(meshX',meshY',reshape(fem_u-disp_u,nx+1,ny+1)');
    title('err-u');
    subplot(2,1,2);
    surf(meshX',meshY',reshape(fem_v-disp_v,nx+1,ny+1)');
    title('err-v');
    % axis equal;
end


