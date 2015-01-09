
function [ ps ] = LagInterpPoly(xs, ys)

n = length(xs);
ps = zeros(n,n);

for i = 1:n
    % the polynomial in the numerator
    pp = poly(xs((1:n)~=i));
    % denominator
    ps(:,i) = pp ./ polyval(pp, xs(i));
end

if (nargin == 2) 
    % if (isrow(ys))
        % ys = ys';
    % end
    % pp = ps * ys;
    % ps = pp;
    
    ps = SumPoly(ps, ys);
end


return
end


