function [fespace] = FESpaceMake(mesh, elem_type, var_dim)

fespace = FESpace();

% fespace.mesh = mesh;

%
fespace.mesh = mesh;

%
fespace.type = elem_type;
fespace.vdim = var_dim;

%
switch elem_type
case 'Q1'
    if mesh.type == 'quad'
        % fespace.nVert = mesh.nVert;
        fespace.nNode = mesh.nVert;
        fespace.nElem = mesh.nElem;
        fespace.nodes = mesh.verts;
        fespace.elems = mesh.conn;
    else
        error('Wrong mesh.type=%s',mesh.type);
    end
otherwise
    error('Unsupport ELEM_TYPE=%s',elem_type);
end

%
fespace.nDofs = fespace.nNode;
fespace.nSize = fespace.nNode * fespace.vdim;



return
end

