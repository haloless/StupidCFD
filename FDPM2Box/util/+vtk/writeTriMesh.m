function [] = writeTriMesh(filename, tri, nodeX,nodeY, datamap)
%writeTriMesh Write triangulated mesh and data. 
% DATAMAP should be a containers.Map.

% useful constants
VTK_TRIANGLE = vtk.vtkCellType('VTK_TRIANGLE');

%
nodeX = nodeX(:);
nodeY = nodeY(:);

fp = fopen(filename, 'w');

fprintf(fp, '# vtk DataFile Version 2.0\n');
fprintf(fp, 'TriData\n');
fprintf(fp, 'ASCII\n');
fprintf(fp, 'DATASET UNSTRUCTURED_GRID\n');

% points
numNodes = length(nodeX);
fprintf(fp, 'POINTS %d float\n', numNodes);
fprintf(fp, '%g %g %g\n', [nodeX,nodeY,zeros(numNodes,1)].');

% cells
numCells = size(tri, 1);
fprintf(fp, 'CELLS %d %d\n', numCells, 4*numCells);
fprintf(fp, '%d %d %d %d\n', [repmat(3, numCells,1), tri-1].'); % NOTE 0-based index
fprintf(fp, 'CELL_TYPES %d\n', numCells);
fprintf(fp, '%d\n', repmat(VTK_TRIANGLE, 1,numCells));


if exist('datamap','var') % write data
    fprintf(fp, 'POINT_DATA %d\n', numNodes);
    names = keys(datamap);
    for i = 1:length(names)
        key = names{i};
        val = datamap(key);
        if size(val,2) == 1
            % scalar field
            fprintf(fp, 'SCALARS %s float 1\n', key);
            fprintf(fp, 'LOOKUP_TABLE default\n');
            fprintf(fp, '%g\n', val.');
        else
            % vector field
            if size(val,2) == 2 % pad z-component by zero
                val = [val, zeros(numNodes,1)];
            end
            fprintf(fp, 'VECTORS %s float\n', key);
            fprintf(fp, '%g %g %g\n', val.');
        end
    end
end

fclose(fp);

return
end



