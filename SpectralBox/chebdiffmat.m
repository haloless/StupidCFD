%%
%%


function [D,x] = chebdiffmat(n)

if n == 0
	D = 0;
	x = 1;
else
	x = chebnode(n);
	
	c = [2; ones(n-1,1); 2] .* (-1).^(0:n)';
	m = repmat(x,1,n+1);
	dm = m - m';
	
	% off-diagonal
	D = (c*(1 ./ c)') ./ (dm+eye(n+1));
	% diagonal
	D = D - diag(sum(D'));
end

return
end


