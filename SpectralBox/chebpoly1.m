
function [cp] = chebpoly(np,xs)

if isrow(xs)
	xs = xs';
end

nx = numel(xs);
% np = np + 1;

cp = zeros(nx,np);

% T0
cp(:,1) = 1;

if np >= 2
	% T1
	cp(:,2) = cos(xs);
	
	if np >= 3
		% n>2, use recurrence relation 
		for i = 3:np
			cp(:,i) = 2*xs.*cp(:,i-1) - cp(:,i-2);
		end
	end
end


return
end

