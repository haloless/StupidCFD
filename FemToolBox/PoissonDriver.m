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

clear;


%% source function

% srctype = 'sinsin'
% srctype = 'quadratic'
% srctype = 'linear'
srctype = 'rz_poly'
switch srctype
case 'sinsin' 
    % f=sin(m*pi*x) * sin(m*pi*y)
    mm = 1; nn = 1;
    % mm = 2; nn = 3;
    fun_fext = @(x,y) (mm^2+nn^2)/2 * sin(mm*pi*x) .* sin(nn*pi*y);
    fun_phi = @(x,y) 1/(2*pi^2) * sin(mm*pi*x) .* sin(nn*pi*y);
    fun_phix = @(x,y) mm/(2*pi) * cos(mm*pi*x) .* sin(nn*pi*y);
    fun_phiy = @(x,y) nn/(2*pi) * sin(mm*pi*x) .* cos(nn*pi*y);
    
case 'quadratic'
    % f = 
    fun_fext = @(x,y) -6.0 * ones(size(x));
    fun_phi  = @(x,y) x.^2 + 2*y.^2 + 1;
    
case 'linear'
    % f = a*x + b*y + c
    aa = 1; bb = 3; cc = 0;
    fun_fext = @(x,y) zeros(size(x));
    fun_phi = @(x,y) cc + aa*x + bb*y;
    fun_phix = @(x,y) aa * ones(size(x));
    fun_phiy = @(x,y) bb * ones(size(y));
    
case 'rz_poly'
    fun_fext = @(r,z) -8.0 * ones(size(r));
    fun_phi = @(r,z) r.^2 + 2*z.^2 + 1;
    
otherwise
    error('unknown srctype: %s', srctype);
end

if strfind(srctype, 'rz_')
    coordtype = 'RZ';
else
    coordtype = 'Cart';
end



%% domain
Lx = 1;
Ly = 1;

% nx = 10;
% ny = 10;
nx = 20;
ny = 20;
% nx = 40;
% ny = 40;
nx1 = nx + 1;
ny1 = ny + 1;

%% build mesh
disp(['Generate mesh']);
mesh = MeshMakeSquare([Lx,Ly], [nx,ny], ...
'xlo',[0,0], 'coord',coordtype);
if 1 % do some transform...
    % mesh.verts = ([1.0,0.2;0.2,1.0] * mesh.verts.').';
    Xs = mesh.verts(:,1);
    Ys = mesh.verts(:,2);
    mesh.verts(:,1) = mesh.verts(:,1) + 0.1 + 0.02*sin(4*pi.*Ys);
    mesh.verts(:,2) = mesh.verts(:,2) + 0.1 + 0.03*sin(2*pi.*Xs);
end

if (1)
    figure;
    MeshPlotEasy(mesh);
    % return
end

%% define FE

% Q1 element, scalar 
fes = FESpaceMake(mesh, 'Q1', 1);
% do not use mesh any more

%% matrix system

% linear form
disp(['Assemble linear']);
fe_b = FELinearForm(fes);
fe_b.setZero();
fe_b.assembleDomain(fun_fext);

% bilinear form
disp(['Assemble bilinear']);
fe_a = FEBilinearForm(fes);
fe_a.assembleDomainLaplacian();


%
numNodes = fes.nNode;
numElems = fes.nElem;
% DOF with [phi]
numDofs = fes.nSize;

nodeX = fes.nodes(:,1);
nodeY = fes.nodes(:,2);
nodes = fes.nodes;
T = fes.elems;

%
Kmat = fe_a.mat;
rhs = fe_b.vec;

%% boundary condition
disp(['Apply BC']);

% analytical solution 
phis = fun_phi(nodeX,nodeY);

disp_ids = unique(mesh.bconn(:));
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
    % return;
end

free_dofs = ones(numDofs,1)';
free_dofs(disp_dofs) = 0;
free_dofs = find(free_dofs);


% set natural BC
rhs(trac_dofs) = trac_vals;

% enforce essential BC
rhs(free_dofs) = rhs(free_dofs) - Kmat(free_dofs,disp_dofs)*disp_vals;

%% 
disp(['Solving equations.']);

sol = zeros(numDofs,1);

sol(free_dofs) = Kmat(free_dofs,free_dofs) \ rhs(free_dofs);
sol(disp_dofs) = disp_vals;


if (1)
    xx = reshape(nodeX, nx1,ny1);
    yy = reshape(nodeY, nx1,ny1);
    
    figure;
    zz = reshape(sol, nx1,ny1);
    % surf(xx,yy,zz,'FaceColor','interp');
    surfc(xx,yy,zz);
    xlabel('x'); ylabel('y'); zlabel('u');
    
    figure;
    phierr = sol-phis;
    surf(xx,yy,reshape(phierr, nx1,ny1));
end




