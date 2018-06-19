function [mesh] = MeshMakeSquare(len,ncell, varargin)

% space dim
ndim = length(len);
if ndim ~= 2
    error('Unsupport NDIM=%d', ndim);
end

% parse varargin
type = 'quad'; % TODO tri quad tetra hex
xlo = zeros(ndim,1);
coord = 'Cart';
for i = 1:2:length(varargin)
    key = varargin{i};
    val = varargin{i+1};
    switch key
    case 'xlo'
        xlo = val;
    case 'coord'
        coord = val;
    otherwise
        error('Unknown %s=%s', key,val);
    end
end

%
x0 = xlo(1);
y0 = xlo(2);
Lx = len(1);
Ly = len(2);
nx = ncell(1);
ny = ncell(2);
nx1 = nx + 1;
ny1 = ny + 1;

%
numElems = nx * ny;
numNodes = nx1 * ny1;

[meshX,meshY] = ndgrid(linspace(x0,x0+Lx,nx1),linspace(y0,y0+Ly,ny1));

nodeX = reshape(meshX,[],1);
nodeY = reshape(meshY,[],1);
nodes = [nodeX, nodeY];

% index function
ind = @(i,j) i + (j-1)*nx1;

% build connectivity
nve = 4; % TODO
T = zeros(numElems,nve);
iElem = 0;
for j = 1:ny
    for i = 1:nx
        iElem = iElem + 1;
        n1 = ind(i,j);
        n2 = ind(i+1,j);
        n3 = ind(i+1,j+1);
        n4 = ind(i,j+1);
        T(iElem,:) = [n1, n2, n3, n4];
    end
end

% build boundary
nbndry = nx*2 + ny*2;
tbndry = [];
tlabel = [];
for i = 1:nx
    tbndry(end+1,:) = [ ind(i,1),ind(i+1,1) ];
    tlabel(end+1,:) = [ 1 ];
end
for j = 1:ny
    tbndry(end+1,:) = [ ind(nx1,j),ind(nx1,j+1) ];
    tlabel(end+1,:) = [ 2 ];
end
for i = nx:-1:1
    tbndry(end+1,:) = [ ind(i+1,ny1),ind(i,ny1) ];
    tlabel(end+1,:) = [ 3 ];
end
for j = ny:-1:1
    tbndry(end+1,:) = [ ind(1,j+1),ind(1,j) ];
    tlabel(end+1,:) = [ 4 ];
end

%
mesh = struct();

%
mesh.nDim = ndim;
mesh.type = type;
mesh.coord = coord;
%
mesh.nElem = numElems;
mesh.nVert = numNodes;
mesh.conn = T;
% mesh.vertX = nodeX;
% mesh.vertY = nodeY;
mesh.verts = nodes;
%
mesh.nBndry = nbndry;
mesh.bconn = tbndry;
mesh.battr = tlabel;

return
end

