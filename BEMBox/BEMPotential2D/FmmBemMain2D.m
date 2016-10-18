
clear variables;

fmm_verbose = 1;

BemMeshGlobals;

FmmTreeGlobals;

%
% maxl: 	    maximum number of elements in a leaf
% levmx:	    maximum number of tree levels
% nexp:	    order of the fast multipole expansions (p)
% ntylr:	    order of the local expansions (= p, in general)
% tolerance:  tolerance for convergence used in the iterative solver
% maxia:	    maximum number of parameters
% ncellmx:    maximum number of cells allowed in the tree
% nleafmx:    maximum number of leaves allowed in the tree
% mxl:	    maximum dimension of Krylov subspace used in the iterative solver
% nwksz:	    size of the space used to store coefficients in preconditioner
% (use default in the code, if value = 0)
%


maxl = 10;
levmx = 10;
nexp = 15;
ntylr = nexp;
tol = 1.0e-8;
maxia = 50;
ncellmx = 5000;
nleafmx = 5000;
mxl = 50;
nwksz = 90000000;

%%
%% create mesh
%%

% input_file = 'input.dat'; 
% input_file = 'input_example_M30.dat'; 
% input_file = 'input_ring_M360q.dat'; 
input_file = 'input_square_M400.dat'; 
BemLoadMesh;


% a         = multipole expansion moments
% b         = local expansion coefficients
% xmax,xmin = maximum and minimum x coordinate
% ymax,ymin = maximum and minimum y coordinate
% ielem     = ielem(i) gives the original element number for i-th element in 
%             the quad-tree structure
% itree     = itree(c) gives the cell location of c-th cell within each tree level
% loct      = elements included in the c-th cell are listed starting at 
%             the loct(c)-th place in the array ielem
% numt      = numt(c) gives the number of elements included in the c-th cell
% ifath     = ifath(c) gives the cell number of the parent cell of the c-th cell
% level     = level l cells start at the level(l)-th cell in the tree
% lowlev    = number of the tree levels


% generate quadtree
FmmGenTree2D;


% multipole moment
a = zeros(nexp+1,ncellmx);      % have zero index
b = zeros(ntylr+1,ncellmx);     % have zero index

% rhs vector
% swap bc
for i = 1:n
    if bc(1,i) == bc_dir
        bc(1,i) = bc_neu;
    elseif bc(1,i) == bc_neu
        bc(1,i) = bc_dir;
    end
end

u = zeros(n,1);
ax = zeros(n,1);
for i = 1:n
    u(i) = bc(2,ielem(i));
end

a = FmmUpward(u,a);
[ax,b] = FmmDownward(u,ax,a,b);

u = -ax;

% rhs = zeros(n,1);
% for i = 1:n
    % rhs(ielem(i)) = u(i);
% end
rhs = u;

for i = 1:n
    if bc(1,i) == bc_dir
        bc(1,i) = bc_neu;
    elseif bc(1,i) == bc_neu
        bc(1,i) = bc_dir;
    end
end

afun = @(xxx) FmmMatVec(xxx,ax,a,b);
[sol,ret,res,iter] = bicgstab(afun,rhs,1.0e-6,2000);
disp(['solver: ret=',int2str(ret),'; res=',num2str(res), '; iter=',int2str(iter)]);

for i = 1:n
    u(ielem(i)) = sol(i);
end


