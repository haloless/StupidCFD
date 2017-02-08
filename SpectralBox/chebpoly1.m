
function [cp] = chebpoly1(np,xs)

if isrow(xs)
	xs = xs';
end

nx = numel(xs);

% [0,...,np]
cp = zeros(nx,np+1);

% T0
cp(:,1) = 1;

if np >= 1
	% T1
	cp(:,2) = xs;
	
	if np >= 2
		% n>1, use recurrence relation 
		for i = 3:np+1
			cp(:,i) = 2*xs.*cp(:,i-1) - cp(:,i-2);
		end
	end
end


return
end

