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

clc;
clear all;

E0 = 500;
nu0 = 0.25;

% plane stress
Dmat = E0 / (1-nu0^2) * [1,nu0,0; nu0,1,0; 0,0,(1-nu0)/2];


L0 = 24;
c0 = 2;
d0 = c0 *2;
I0 = d0^3 / 12;
% end load
P0 = 1;

elem_type = 'Q4';
elem_numNodes = 4;
elem_numGauss = 4;

refine = 2;
nx = refine * 24;
ny = refine * 4;

% TODO move the following operations to GenerateMesh
% generate mesh
numElems = nx * ny;
numNodes = (nx+1) * (ny+1);
[meshX,meshY] = ndgrid(linspace(0,L0,nx+1),linspace(-c0,c0,ny+1));
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
timo_usol = @(x,y) -P0/(6*E0*I0)*y .* ((6*L0-3*x).*x + (2+nu0)*(y.^2-c0^2));
timo_vsol = @(x,y) P0/(6*E0*I0) .* (3*nu0*y.^2.*(L0-x) + c0^2*(4+5*nu0)*x + (3*L0-x).*x.^2);
timo_u = timo_usol(nodeX,nodeY);
timo_v = timo_vsol(nodeX,nodeY);



% global stiffness matrix
Kmat = sparse(numNodes*2,numNodes*2);
rhs = zeros(numNodes*2,1);

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
% disp_ids = [1:nx+1:numNodes, nx+1:nx+1:numNodes]';
disp_ids = [1:nx+1:numNodes]';
disp_dofs = zipcol([disp_ids*2-1, disp_ids*2]);
disp_vals = zipcol([timo_u(disp_ids),timo_v(disp_ids)]);

trac_ids = [nx+1+(ny/2)*(nx+1)];
trac_dofs = zipcol([trac_ids*2-1, trac_ids*2]);
trac_vals = [0; P0];
if (1)
    figure;
    plot(nodeX(disp_ids),nodeY(disp_ids),'x', nodeX(trac_ids),nodeY(trac_ids),'o'); 
    legend('disp','trac');
    axis equal;
end

% enforce trac. BC
rhs(trac_dofs) = trac_vals;
% enforce disp. BC
spI = speye(numNodes*2,numNodes*2);
Kmat(disp_dofs,:) = spI(disp_dofs,:);
rhs(disp_dofs) = disp_vals;

disp(['Solving equations.']);
sol = Kmat \ rhs;
sol = reshape(sol,2,[]);
fem_u = sol(1,:)';
fem_v = sol(2,:)';

if (1)
    figure;
    plot(nodeX+fem_u,nodeY+fem_v,'o', nodeX+timo_u,nodeY+timo_v,'x'); 
    legend('FEM','timo'); axis equal;
    
    figure;
    subplot(2,1,1);
    surf(meshX',meshY',reshape(fem_u-timo_u,nx+1,ny+1)');
    title('err-u');
    subplot(2,1,2);
    surf(meshX',meshY',reshape(fem_v-timo_v,nx+1,ny+1)');
    title('err-v');
    % axis equal;
end


