
function [ coef ] = LagInterpCoef(xint, xs)

coef = zeros(size(xs));

n = length(xs);

for i = 1:n
    % js = [1:i-1, i+1:n];
    js = ((1:n) ~= i);
    coef(i) = prod(xint-xs(js)) / prod(xs(i)-xs(js));
end

return
end


