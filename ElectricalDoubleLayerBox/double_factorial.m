function [val] = double_factorial(n)
% The n!! double factorial
% n!! = prod(1, 3, 5, ..., n)
%

if n == 0
	val = 1;
else
	ns = double(1:2:n);
	val = prod(ns);
end

return
end

