%Load GMSH file

[groups,nodes] = gmsh.loadMeshFile('testBarNecking.msh');
numNodes = size(nodes,1);
nodeX = nodes(:,1);
nodeY = nodes(:,2);

loaded_nodes = [];
loaded_force = [];
fixed_nodes_x = zeros(numNodes,1);
fixed_disp_x = zeros(numNodes,1);
fixed_nodes_y = zeros(numNodes,1);
fixed_disp_y = zeros(numNodes,1);

ids = groups('left').nodeIds;
fixed_nodes_x(ids) = 1;
fixed_disp_x(ids) = 0;

ids = groups('bottom').nodeIds;
fixed_nodes_y(ids) = 1;
fixed_disp_y(ids) = 0;

ids = groups('top').nodeIds;
fixed_nodes_x(ids) = 1;
fixed_disp_x(ids) = 0;
fixed_nodes_y(ids) = 1;
fixed_disp_y(ids) = Ldisp;

%
fixed_nodes_x = find(fixed_nodes_x==1);
fixed_disp_x = fixed_disp_x(fixed_nodes_x);
fixed_nodes_y = find(fixed_nodes_y==1);
fixed_disp_y = fixed_disp_y(fixed_nodes_y);

% setup Dirichlet BC
tmp = fdpmNodeDof(fixed_nodes_y,'y');
dispBCDofs = [dispBCDofs; tmp];
dispBCVals = [dispBCVals; fixed_disp_y];
tmp = fdpmNodeDof(fixed_nodes_x,'x');
dispBCDofs = [dispBCDofs; tmp];
dispBCVals = [dispBCVals; fixed_disp_x];

% setup Loading
tracBCDofs = fdpmNodeDof(loaded_nodes,'xy');
tracBCVals = reshape(loaded_force, [],1);


%
edgePair = [];
edgePair = [edgePair; groups('left').conn];
edgePair = [edgePair; groups('bottom').conn];
edgePair = [edgePair; groups('right').conn];
edgePair = [edgePair; groups('top').conn];
numEdges = size(edgePair,1);


