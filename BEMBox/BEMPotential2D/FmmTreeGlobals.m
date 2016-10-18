
%%
%% maxl: 	    maximum number of elements in a leaf
%% levmx:	    maximum number of tree levels
%% nexp:	    order of the fast multipole expansions (p)
%% ntylr:	    order of the local expansions (= p, in general)
%% tolerance:  tolerance for convergence used in the iterative solver
%% maxia:	    maximum number of parameters
%% ncellmx:    maximum number of cells allowed in the tree
%% nleafmx:    maximum number of leaves allowed in the tree
%% mxl:	    maximum dimension of Krylov subspace used in the iterative solver
%% nwksz:	    size of the space used to store coefficients in preconditioner
%% (use default in the code, if value = 0)
%%

global maxl;
global levmx;
global nexp ntylr;
global ncellmx nleafmx;

% ielem     = ielem(i) gives the original element number for i-th element in 
%             the quad-tree structure
% itree     = itree(c) gives the cell location of c-th cell within each tree level
% loct      = elements included in the c-th cell are listed starting at 
%             the loct(c)-th place in the array ielem
% numt      = numt(c) gives the number of elements included in the c-th cell
% ifath     = ifath(c) gives the cell number of the parent cell of the c-th cell
% level     = level l cells start at the level(l)-th cell in the tree
% lowlev    = number of the tree levels
global ielem;
global itree;
global level;
global loct;
global numt;
global ifath;
% leaf=1, otherwise=0
global leafflag;
% the finest level
global lowlev;

% for easy access
global levelndiv;
global leveldx leveldy;







