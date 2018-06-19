function [n] = voigt_norm(v)

if numel(v) ~= 6
	error('input vector size is not 6')
end

va = v(1:3);
vb = v(4:6);
n = sqrt(sum(va.^2 + 2*(vb.^2)));


return
end

