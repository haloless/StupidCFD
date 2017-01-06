% may also use matlab's PINV
function [minv] = SafeInvert(mat)

% tol = epsilon * max(m,n) * max(Sigma)

esmall = 1.0e-6;

test = det(mat);

if abs(test) >= esmall
	minv = inv(mat);
else
	% near-sigular, use SVD to obtain pseudo-inverse
	[u,s,v] = svd(mat);
	sinv = s';
	n = size(mat,1);
	for i = 1:n
		if abs(s(i,i)) > esmall
			sinv(i,i) = 1.0 / s(i,i);
		else
			sinv(i,i) = 0.0;
		end
	end
	
	minv = v * sinv * u';
end

return
end

