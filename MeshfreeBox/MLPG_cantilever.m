% ## Copyright (C) 2013 admin_2
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

% ## MLPG_cantilever

% ## Author: admin_2 <admin_2@KOSHIZUKA>
% ## Created: 2013-07-22

clc;
clear all;

D = 4;
L = 24;
E = 1000;
nu = 0.25;
P = 1;
I = D^3 / 12;
% plain stress material
Emat = E/(1-nu^2) * [1,nu,0; nu,1,0; 0,0,(1-nu)/2];
Emat = sparse(Emat);
% analytical solution
timo_u = @(x,y) -P/(6*E*I) * y .* ((6*L-3*x).*x + (2+nu)*(y.^2-D^2/4));
timo_v = @(x,y) P/(6*E*I) * (3*nu * y.^2 .* (L-x) + (4+5*nu)*D^2*x/4 + (3*L-x).*x.^2);
timo_sxx = @(x,y) -P/I * (L-x) .* y;
timo_syy = @(x,y) zeros(size(x));
timo_sxy = @(x,y) P/(2*I) * (D^2/4 - y.^2);

refine = 2;
nx = refine*12 + 1;
ny = refine*2 + 1;
N = nx * ny;
hx = L / (nx-1);
hy = D / (ny-1);
hh = max([hx,hy]);
dilation = 4.0;
re = hh * dilation;

[Xs,Ys] = ndgrid((0:nx-1)*hx,(0:ny-1)*hy-D/2);
% figure; plot(Xs',Ys',"x"); axis equal; axis([0, L, -D/2, D/2]);
% figure; plot((Xs+timo_u(Xs,Ys))',(Ys+timo_v(Xs,Ys))',"x"); axis equal; %axis([0, L, -D/2, D/2]);
xs = reshape(Xs,[],1);
ys = reshape(Ys,[],1);

trac_nodes = [nx:nx:N, 2:nx-1, (ny-1)*nx+2:(ny-1)*nx+nx-1];
% trac_nodes = [nx:nx:N];
% trac_nodes = [];
disp_nodes = 1:nx:N;
% disp_nodes = [1:nx:N, nx:nx:N];
if (0)
	figure; 
	plot(Xs(disp_nodes)',Ys(disp_nodes)','x',Xs(trac_nodes)',Ys(trac_nodes)','o'); 
	axis equal;
	legend('displacement','traction');
end
trac_dofs = reshape([trac_nodes*3-2; trac_nodes*3-1; trac_nodes*3],1,[]);
trac_vals = reshape([timo_sxx(Xs(trac_nodes),Ys(trac_nodes)); timo_syy(Xs(trac_nodes),Ys(trac_nodes)); timo_sxy(Xs(trac_nodes),Ys(trac_nodes))],1,[]);
disp_dofs = reshape([disp_nodes*2-1; disp_nodes*2],1,[]);
disp_vals = reshape([timo_u(Xs(disp_nodes),Ys(disp_nodes)); timo_v(Xs(disp_nodes),Ys(disp_nodes))],1,[]);
disp_dofs = disp_dofs + N*3; % global offset between sigma and u/v

% build MLS
MLSNeigh = sparse(N,N);
MLSShape = sparse(N,N);
MLSGradx = sparse(N,N);
MLSGrady = sparse(N,N);
for i = 1:N
	xc = xs(i);
	yc = ys(i);
	
	% connectivity
	neigh = MLS_neigh(xc,yc,xs,ys,re);
	MLSNeigh(i,:) = neigh;
	
	% shape func.
	[phi,phi_x,phi_y] = MLS_shape(xc,yc,xs,ys,re,neigh);
	MLSShape(i,:) = phi;
	MLSGradx(i,:) = phi_x;
	MLSGrady(i,:) = phi_y;
end
disp(['MLS parts constructed.']);

% work in true nodal value space
if (1)
	MT = MLSShape';
	MLSShape = speye(N,N);
	MLSGradx = (MT \ MLSGradx')';
	MLSGrady = (MT \ MLSGrady')';
	disp(['MLS transformed to true nodal value space.']);
end

% construct stiffness and constitution matrix
Ks = sparse(N*2,N*3);
T = sparse(N*3,N*2);
for i = 1:N
	phix = MLSGradx(i,:);
	phiy = MLSGrady(i,:);
	
	Ks(i*2-1,:) = reshape([phix;zeros(1,N);phiy],1,[]);
	Ks(i*2,:) = reshape([zeros(1,N);phiy;phix],1,[]);
	
	Gi = [reshape([phix;zeros(1,N)],1,[]); reshape([zeros(1,N);phiy],1,[]); reshape([phiy;phix],1,[])];
	T(i*3-2:i*3,:) = sparse(Emat * Gi);
	% T(i*3-2,:) = E/(1-nu^2) * reshape([phix;nu*phiy],1,[]);
	% T(i*3-1,:) = E/(1-nu^2) * reshape([nu*phix;phiy],1,[]);
	% T(i*3,:) = E/(1-nu^2) * (1-nu)/2 * reshape([phiy;phix],1,[]);
end
disp(['Ks & T constructed.']);

% global matrix
K = [Ks, sparse(2*N,2*N); speye(3*N,3*N), -T];
% assume zero body force
rhs = zeros(3*N+2*N,1);

% enforce BC
K(trac_dofs,:) = speye(5*N)(trac_dofs,:);
rhs(trac_dofs) = trac_vals;

K(disp_dofs,:) = speye(5*N)(disp_dofs,:);
rhs(disp_dofs) = disp_vals;

% solution
sol = K \ rhs;
stress = reshape(sol(1:3*N),3,[]);
disp = reshape(sol(3*N+1:5*N),2,[]);
u = disp(1,:)'; v = disp(2,:)';

% true solution
timou = timo_u(xs,ys);
timov = timo_v(xs,ys);
timodisp = reshape([timou';timov'],[],1);
timos11 = timo_sxx(xs,ys);
timos22 = timo_syy(xs,ys);
timos12 = timo_sxy(xs,ys);
timostress = reshape([timos11';timos22';timos12'],[],1);

%
figure;
%
sigma = T * timodisp;
sigma = reshape(sigma,3,[]);
subplot(2,3,1); plot(sigma(1,:),"x", timos11,"-"); title("s11-timo");
subplot(2,3,2); plot(sigma(2,:),"x", timos22,"-"); title("s22-timo");
subplot(2,3,3); plot(sigma(3,:),"x", timos12,"-"); title("s12-timo");
%
sigma = T * disp(:);
sigma = reshape(sigma,3,[]);
sigma = stress;
subplot(2,3,4); plot(sigma(1,:),"x", timos11,"-"); title("s11-numer");
subplot(2,3,5); plot(sigma(2,:),"x", timos22,"-"); title("s22-numer");
subplot(2,3,6); plot(sigma(3,:),"x", timos12,"-"); title("s12-numer");

if (0)
	fb = reshape(Ks*stress(:),2,[]);
	figure; 
	plot(fb(:),"x", Ks*timostress,"-"); title("residual force");
end

figure;
plot(xs+u,ys+v,'x', xs+timou,ys+timov,'o');
axis equal;
legend("numer.", "Timo.");

verr = mean(abs(v-timov))


