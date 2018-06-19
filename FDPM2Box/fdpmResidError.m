function [err] = fdpmResidError(resid, measure)
%fdpmResidError
% Calculates the residual error
% Default is L2-error
%

% default error measure
if nargin == 1
	measure = 'L2';
end

% ensure vector
resid = resid(:);


switch measure
	case 'L2'
		err = sqrt(mean(resid.^2));
	case 'L1'
		err = mean(abs(resid));
	case 'Linf'
		err = max(abs(resid));
	otherwise
		error(['Unknown error measure = ', measure]);
end



return
end


