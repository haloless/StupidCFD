function [vec] = voigtMat2Vec(mat,op)
%voigtMat2Vec
% Voigt notation, encode

if numel(mat) ~= 9
	error('Invalid input matrix');
end

switch op
	case 'stress'
		vec = [mat(1); mat(5); mat(9);   mat(2);   mat(6);   mat(7)];
	case 'strain'
		vec = [mat(1); mat(5); mat(9); 2*mat(2); 2*mat(6); 2*mat(7)];
	otherwise
		error(['Unknown voigt type=',op]);
end


return
end

