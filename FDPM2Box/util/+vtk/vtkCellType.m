function [vtype] = vtkCellType(key)
%vtkCellType Load constants for VTK cell type.


vtype = -1;

switch key
case 'VTK_VERTEX'
    vtype = 1;
case 'VTK_POLY_VERTEX'
    vtype = 2;
case 'VTK_LINE'
    vtype = 3;
case 'VTK_POLY_LINE'
    vtype = 4;
case 'VTK_TRIANGLE'
    vtype = 5;
case 'VTK_TRIANGLE_STRIP'
    vtype = 6;
case 'VTK_POLYGON'
    vtype = 7;
case 'VTK_PIXEL'
    vtype = 8;
case 'VTK_QUAD'
    vtype = 9;
case 'VTK_TETRA'
    vtype = 10;
case 'VTK_VOXEL'
    vtype = 11;
case 'VTK_HEXAHEDRON'
    vtype = 12;
case 'VTK_WEDGE'
    vtype = 13;
case 'VTK_PYRAMID'
    vtype = 14;
case 'VTK_QUADRATIC_EDGE'
    vtype = 21;
case 'VTK_QUADRATIC_TRIANGLE'
    vtype = 22;
case 'VTK_QUADRATIC_QUAD'
    vtype = 23;
case 'VTK_QUADRATIC_TETRA'
    vtype = 24;
case 'VTK_QUADRATIC_HEXAHEDRON'
    vtype = 25;
otherwise
    error('Unknown cell type=%s', key);
end



return
end


