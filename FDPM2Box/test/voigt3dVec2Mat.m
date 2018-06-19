function [mat] = voigtVec2Mat(vec,op)
%voigtVec2Mat
% Voigt notation, decode

if numel(vec) ~= 6
	error('Invalid input vector');
end

switch op
	case 'stress'
		mat = [vec(1),vec(4),vec(6); vec(4),vec(2),vec(5); vec(6),vec(5),vec(3)];
	case 'strain'
		mat = [vec(1),vec(4)/2,vec(6)/2; vec(4)/2,vec(2),vec(5)/2; vec(6)/2,vec(5)/2,vec(3)];
	otherwise
		error(['Unknown voigt type=',op]);
end


return
end

