function [v,d,repeat] = mathEigen(a)
%mathEigen: similar to [v,d] = eig(a)
% Only operates on 2x2
% Checks repeated eigenvalues

[v,d] = eig(a(1:2,1:2));

if nargout >= 3
	% check repeated eigenvalues
	% tol_abs = 1.0e-9;
	% tol_rel = 1.0e-6;
	tol = 1.0e-6;
	
	repeat = 0;
	d1 = d(1,1);
	d2 = d(2,2);
	diff = abs(d1-d2);
	dmax = max(abs(d1), abs(d2));
	if dmax > 0
		diff = diff / dmax;
	end
	if diff < tol
		repeat = 1;
	end
end



return
end


