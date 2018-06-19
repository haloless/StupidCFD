function [y] = fdpmLogm(x)
%fdpmLogm
% special y = logm(x)
% for x is a plane-strain or axisymmetric 3x3 matrix
% [ x11, x12, 0   ]
% [ x12, x22, 0   ]
% [ 0,   0,   x33 ]
%

y = zeros(3);

if 0
	% matlab logm is not very fast...
	y(1:2,1:2) = logm(x(1:2,1:2));
end

if 1
	% eig is faster for small 2x2 submatrix
	[v,d] = eig(x(1:2,1:2));
	y(1:2,1:2) = v * diag(log(diag(d))) * v.';
end



% the 33 component
y(3,3) = log(x(3,3));


return
end


