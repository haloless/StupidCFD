function [vproj] = mathEigenProj(v,d,repeat)
%mathEigenProj: eigenprojection from [V,D]
%

if repeat
	warning('Repeated eigenvalues, eigenprojection incorrect!');
end

vproj = zeros(4,2);
for dir = 1:2
	vv = v(:,dir) * v(:,dir)';
	vproj(1:3,dir) = [vv(1,1),vv(2,2),vv(1,2)];
end


return
end

