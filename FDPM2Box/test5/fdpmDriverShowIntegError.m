%fdpmDriverShowIntegError: Check integral residual of shape derivative
% OUTPUT
% - v_sum
% 

%% NOTE on axisymmetric RZ coord
%
% $$\nabla^2 \phi + f = 0$$
%
% $$\frac{1}{r} \frac{\partial}{\partial r} \left( r \frac{\partial \phi}{\partial r}\right) + \frac{\partial^2 \phi}{\partial z^2}$$
%
% In RZ coord, a linear function $\phi = ar + bz + c$.
% This gives $\phi_r = a$, $\phi_z = b$,
% hence $f = -a/r$.
% 
% The weak form is 
%
% $$\int f v r dr dz + \int v \mathbf{q} \cdot \mathbf{n} r ds = \int \nabla u \cdot \nabla v r dr dz$$
% 
% And the integral consistency requires 
%
% $$
% \int \nabla v rdrdz + \int \pmatrix{1 \cr 0} v drdz
% = \int v \mathbf{n} rds
% $$
%
% This finally becomes
%
% $$\int \nabla (vr) drdz = \int v \mathbf{n} rds$$
%
%



%
v_sum = zeros(numNodes,2);

% 
for i = 1:numNodes
    ineigh = conn(i).neigh2;
    ivol = nodeVol(i);
    
    v_sum(ineigh,1) = v_sum(ineigh,1) + ivol.*conn(i).dNX(:);
    v_sum(ineigh,2) = v_sum(ineigh,2) + ivol.*conn(i).dNY(:);
end

% special hoop term for RZ coord
if exist('prob_type','var') && prob_type==2
    for i = 1:numNodes
        ineigh = conn(i).neigh2;
        iarea = nodeArea(i);
        
        v_sum(i,1) = v_sum(i,1) + 1 * iarea;
    end
end

figure; 
plot(nodeX, nodeY, 'k.');
axis equal;
hold on;
quiver(nodeX(:),nodeY(:),v_sum(:,1),v_sum(:,2));
hold off;
