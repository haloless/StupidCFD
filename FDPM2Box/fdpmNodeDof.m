function [idof] = fdpmNodeDof(inode,op)
%fdpmNodeDof
% Node index -> DOF index
% Input
% OP: 0->xy; 1->x; 2->y;
% Output
% column vector idof

% ensure a row vector
inode = inode(:).';

switch op
	case 'xy'
		idof = reshape([inode*2-1; inode*2], [],1);
	case 'x'
		idof = reshape(inode*2-1, [],1);
	case 'y'
		idof = reshape(inode*2,   [],1);
	otherwise
		disp('Unknown op = ');
		disp(op);
		error('Unknown op');
end



return
end


