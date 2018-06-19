function [matl] = materialNeoHookeanPlaneStrain(E0,nu0)
% materialNeoHookeanPlaneStrain

%
I = eye(2);

% Lame's constants
[lambda0,mu0] = materialLameConst(E0,nu0);


matl = struct();
matl.sigma = @func_sigma;
matl.D = @func_D;
matl.S = @func_S;

%
% the following functions are nested to give enclosures for constitutive equations.
% in principle, J = det(F)
%

	% nested function 
	% cauchy stress
	function [sigma] = func_sigma(Fi,Ji)
		sigma = mu0/Ji*(Fi*Fi'-I) + lambda0*log(Ji)/Ji*I;
	return
	end
	
	% nested function
	% tangent matrix
	function [Dmat] = func_D(Fi,Ji)
		% use modified constants to simplify 
		lambdai = lambda0 / Ji;
		mui = (mu0-lambda0*log(Ji)) / Ji;
		% this is the constitutive matrix
		Dmat = [lambdai+2*mui, lambdai, 0; ...
		lambdai, lambdai+2*mui, 0; ...
		0, 0, mui];
	return
	end
	
	% nested function
	% 2nd PK stress
	function [Si] = func_S(Fi,Ji)
		invCi = inv(Fi' * Fi);
		Si = mu0*(I-invCi) + lambda0*log(Ji)*invCi;
	return
	end

return
end


